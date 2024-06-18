import logging
from collections import OrderedDict
from itertools import groupby, islice
from operator import itemgetter
from typing import Any, Dict, Tuple, Union
from uuid import UUID

import numpy as np
import pandas as pd
from django.core.cache import cache

from ..models import Analysis
from ..utils import clear_data, data_to_edges

logger = logging.getLogger(__name__)


def rename_tf(name: Tuple[str, str, str], gene_name_symbol: str) -> str:
    return f'{gene_name_symbol} "{name[1]}" {name[2]}'


def get_summary(uid: Union[UUID, str], size_limit: int = 50) -> Dict[str, Any]:
    try:
        cached_data = cache.get_many([f"{uid}/tabular_output", f"{uid}/analysis_ids"])
        df, ids = itemgetter(f"{uid}/tabular_output", f"{uid}/analysis_ids")(cached_data)
    except KeyError as e:
        raise ValueError('data not found') from e

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
        ((a.pk, a.name, a.tf.gene_name_symbol) for a in analyses.iterator(chunk_size=2000)),
        columns=['pk', 'name', 'symbol']
    ).set_index('pk')

    df = (df
          .pipe(data_to_edges)
          .apply(lambda x: x.value_counts())
          .fillna(0)
          .astype(np.int_))

    tf_order = df.groupby(axis=1, level=0).sum().sum(axis=0).sort_values(ascending=False)
    tf_total = tf_order.groupby(by=list(map(itemgetter(0), tf_order.index))).sum()
    tf_reorder = sorted(tf_order.index, key=lambda i: (tf_total[i[0]], tf_order.at[i]), reverse=True)

    chart = (df.reindex(columns=tf_reorder, level=0))
    chart.columns = chart.columns.to_flat_index()

    tf_names = {
        (t, a): rename_tf(t, analysis_data.at[a, 'symbol']) + ',' + analysis_data.at[a, 'name'] + (
            '\n' + ids[(t, a)]['name'] if ids[(t, a)]['version'] else '')
        for t, a in chart.columns}

    chart = (chart
             .rename(columns=tf_names)
             .to_dict(into=OrderedDict))

    result = {
        'chart': chart
    }

    if errors:
        result['errors'] = errors

    return result
