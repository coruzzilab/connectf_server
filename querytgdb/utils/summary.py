import re
import time
from collections import OrderedDict
from typing import Any, Dict
import logging

import numpy as np
import pandas as pd

from ..models import Analysis
from ..utils import data_to_edges

logger = logging.getLogger(__name__)


def rename_tf(name: str, gene_id: str, gene_name_symbol: str) -> str:
    return re.sub(r'^' + re.escape(gene_id), gene_name_symbol, name, flags=re.I)


def get_summary(df: pd.DataFrame) -> Dict[str, Any]:
    df = df.loc[:, (slice(None), slice(None), ['EDGE', 'Log2FC'])]

    analyses = Analysis.objects.filter(
        pk__in=df.columns.get_level_values(1)
    ).distinct().prefetch_related('tf')

    analysis_data = pd.DataFrame(
        ((a.pk, a.name, a.tf.gene_id, a.tf.gene_name_symbol) for a in analyses.iterator()),
        columns=['pk', 'name', 'gene_id', 'symbol']
    ).set_index('pk')

    col_names = dict(analysis_data[['name']].itertuples(name=None))
    tf_names = {t: rename_tf(t, analysis_data.at[a, 'gene_id'], analysis_data.at[a, 'symbol']) for t, a, c in
                df.columns}

    df = data_to_edges(df)

    df = df.apply(lambda x: x.value_counts()).fillna(0).astype(np.int_)

    chart = (df.reindex(columns=df.sum(axis=1, level=0).sum(axis=0).sort_values(ascending=False).index, level=0)
             .rename(columns=col_names, level=1)
             .rename(columns=tf_names, level=0)
             .to_dict(into=OrderedDict))

    chart = OrderedDict((','.join(key), value) for key, value in chart.items())

    return {
        'chart': chart
    }
