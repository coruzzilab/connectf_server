from collections import OrderedDict

import numpy as np
import pandas as pd

from ..models import Analysis
from ..utils import data_to_edges


def get_summary(cache_path):
    df = pd.read_pickle(cache_path)

    df = df.loc[:, (slice(None), slice(None), ['EDGE', 'Log2FC'])]

    analyses = Analysis.objects.filter(
        pk__in=df.columns.get_level_values(1)
    ).prefetch_related('analysisdata_set', 'analysisdata_set__key')

    col_names = {a.pk: a.name for a in analyses}

    df = data_to_edges(df, analyses)

    df = df.apply(lambda x: x.value_counts()).fillna(0).astype(np.int_)

    chart = (df.reindex_axis(df.sum(axis=1, level=0).sum(axis=0).sort_values(ascending=False).index, axis=1, level=0)
             .rename(columns=col_names, level=1)
             .to_dict(into=OrderedDict))

    chart = OrderedDict((','.join(key), value) for key, value in chart.items())

    return {
        'chart': chart
    }
