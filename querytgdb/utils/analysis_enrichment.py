from collections import OrderedDict
from itertools import combinations
from typing import Dict, Tuple

import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import fdrcorrection

from ..models import Analysis
from ..utils import split_name
from ..utils.parser import ANNOTATIONS


class AnalysisEnrichmentError(ValueError):
    pass


def split_col_name(col_name: Tuple[str, int]) -> Tuple[str, str, int]:
    name, analysis_id = col_name

    return split_name(name) + (analysis_id,)


def make_col_tuple(col_name: Tuple[str, int], analysis: Analysis) -> Tuple[str, str, int]:
    return (analysis.tf.gene_name_symbol, *split_col_name(col_name)[1:])


def analysis_enrichment(cache_path) -> Dict:
    df = pd.read_pickle(cache_path)

    df = df.loc[:, (slice(None), slice(None), ['EDGE', 'Log2FC'])]
    df.columns = df.columns.droplevel(2)

    if df.shape[1] < 2:
        raise AnalysisEnrichmentError('Analysis enrichment requires more than 1 queried analysis')

    columns = []
    data = []
    info = []

    background = ANNOTATIONS.shape[0]

    analysis_ids = df.columns.get_level_values(1)

    analyses = Analysis.objects.filter(pk__in=analysis_ids).prefetch_related('analysisdata_set', 'tf')

    for col_name in df.columns:
        analysis = analyses.get(pk=col_name[1])
        d = OrderedDict(Count=df[col_name].count())
        d.update(analysis.meta_dict)

        info.append((make_col_tuple(col_name, analysis), d))

    for (name1, col1), (name2, col2) in combinations(((name, col.index[col.notna()])
                                                      for name, col in df.iteritems()), 2):
        columns.append((
            make_col_tuple(name1, analyses.get(pk=name1[1])),
            make_col_tuple(name2, analyses.get(pk=name2[1]))
        ))

        common = col1.intersection(col2).sort_values()

        c = (
            (len(common), len(col1.difference(col2))),
            (len(col2.difference(col1)), background - len(col1.union(col2)))
        )

        data.append({
            'greater': fisher_exact(c, 'greater')[1],
            'less': fisher_exact(c, 'less')[1],
            'genes': common
        })

    p_values = (pd.DataFrame(data)[['less', 'greater']]
                .apply(lambda x: fdrcorrection(x)[1])
                .rename(columns=lambda x: x + '_adj'))

    for d, adj_p in zip(data, p_values.itertuples(index=False)):
        d.update(adj_p._asdict())

    return {
        'columns': columns,
        'data': data,
        'info': info
    }
