from itertools import combinations, islice
from typing import Dict

import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import fdrcorrection

from ..models import Analysis
from ..utils.parser import ANNOTATIONS


def analysis_enrichment(cache_path) -> Dict:
    df = pd.read_pickle(cache_path)

    df = df.loc[:, (slice(None), slice(None), slice(None), 'EDGE')]
    df.columns = df.columns.droplevel(3)

    if df.shape[1] < 2:
        raise ValueError('Analysis enrichment requires more than 1 queried analysis')

    columns = []
    data = []
    info = []

    background = ANNOTATIONS.shape[0]

    exp_ids, analysis_ids = islice(zip(*df.columns), 1, None)

    analyses = Analysis.objects.filter(name__in=analysis_ids, experiment__name__in=exp_ids).prefetch_related(
        'experiment', 'experiment__experimentdata_set', 'analysisdata_set')

    for name in df.columns:
        analysis = analyses.get(name=name[2], experiment__name=name[1])

        i = dict(analysis.analysisdata_set.values_list('key', 'value'))
        i.update(analysis.experiment.experimentdata_set.values_list('key', 'value'))

        info.append((name, i))

    for (name1, col1), (name2, col2) in combinations(((name, col.index[col.notna()])
                                                      for name, col in df.iteritems()), 2):
        columns.append((name1, name2))

        c = (
            (len(col1.intersection(col2)), len(col1.difference(col2))),
            (len(col2.difference(col1)), background - len(col1.union(col2)))
        )

        data.append({
            'greater': fisher_exact(c, 'greater')[1],
            'less': fisher_exact(c, 'less')[1]
        })

    p_values = pd.DataFrame(data).apply(lambda x: fdrcorrection(x)[1]).rename(columns=lambda x: x + '_adj')

    for d, adj_p in zip(data, p_values.itertuples(index=False)):
        d.update(adj_p._asdict())

    return {
        'columns': columns,
        'data': data,
        'info': info
    }
