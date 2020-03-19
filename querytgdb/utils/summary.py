import logging
from collections import OrderedDict
from itertools import groupby, islice
from operator import itemgetter
from typing import Any, Dict, Tuple

import numpy as np
import pandas as pd

from ..models import Analysis
from ..utils import clear_data, data_to_edges

logger = logging.getLogger(__name__)


def rename_tf(name: Tuple[str, str, str], gene_name_symbol: str) -> str:
    return f'{gene_name_symbol} "{name[1]}" {name[2]}'


def get_summary(df: pd.DataFrame, size_limit: int = 50) -> Dict[str, Any]:
    df = clear_data(df, False)

    errors = []

    top_tfs = list(islice(map(itemgetter(0), groupby(df.columns.get_level_values(0))), 0, size_limit))

    if len(top_tfs) == size_limit:
        df = df.loc[:, (top_tfs, slice(None))]
        errors.append(f'Only showing top {size_limit} TFs based on edge count!')

    analyses = Analysis.objects.filter(
        pk__in=df.columns.get_level_values(1)
    ).distinct().prefetch_related('tf')

    analysis_data = pd.DataFrame(
        ((a.pk, a.name, a.tf.gene_name_symbol) for a in analyses.iterator()),
        columns=['pk', 'name', 'symbol']
    ).set_index('pk')

    col_names = dict(analysis_data[['name']].itertuples(name=None))
    tf_names = {t: rename_tf(t, analysis_data.at[a, 'symbol']) for t, a, c in
                df.columns}

    df = data_to_edges(df)

    df = df.apply(lambda x: x.value_counts()).fillna(0).astype(np.int_)

    chart = (df.reindex(columns=df.sum(axis=1, level=0).sum(axis=0).sort_values(ascending=False).index, level=0)
             .rename(columns=col_names, level=1)
             .rename(columns=tf_names, level=0)
             .to_dict(into=OrderedDict))

    chart = OrderedDict((','.join(key), value) for key, value in chart.items())

    result = {
        'chart': chart
    }

    if errors:
        result['errors'] = errors

    return result
