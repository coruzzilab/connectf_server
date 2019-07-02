import csv
import warnings
from collections import OrderedDict
from io import StringIO
from itertools import combinations
from operator import itemgetter
from typing import Dict, Optional, Tuple, Union
from uuid import UUID

import pandas as pd
from django.core.cache import cache
from django.http import HttpResponse
from scipy.special import comb
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import fdrcorrection

from querytgdb.utils import annotations
from ..models import Analysis
from ..utils import clear_data, split_name


class AnalysisEnrichmentError(ValueError):
    pass


class AnalysisEnrichmentWarning(AnalysisEnrichmentError, UserWarning):
    pass


def split_col_name(col_name: Tuple[str, int]) -> Tuple[str, str, int]:
    name, analysis_id = col_name

    return split_name(name) + (analysis_id,)


def make_col_tuple(col_name: Tuple[str, int], analysis: Analysis) -> Tuple[str, str, int]:
    return (analysis.tf.gene_name_symbol, *split_col_name(col_name)[1:])


def analysis_enrichment(uid: Union[UUID, str], size_limit: int = 100, raise_warning: bool = False) -> Dict:
    df = cache.get(f'{uid}/tabular_output')

    if df is None:
        raise AnalysisEnrichmentError("Please make a new query")

    df = clear_data(df)

    if df.shape[1] < 2:
        raise AnalysisEnrichmentError('Analysis enrichment requires more than 1 queried analysis')

    if comb(df.shape[1], 2) > size_limit:
        e = AnalysisEnrichmentWarning('Data size too large.')
        if raise_warning:
            raise e
        else:
            warnings.warn(e)

    columns = []
    data = []
    info = []

    background = annotations().shape[0]

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


def analysis_enrichment_csv(uid: Union[str, UUID],
                            buffer: Optional[Union[StringIO, HttpResponse]] = None,
                            size_limit: int = 100,
                            raise_warning: bool = False):
    if buffer is None:
        buffer = StringIO()

    enrichment = cache.get(f'{uid}/analysis_enrichment')

    if enrichment is None:
        enrichment = analysis_enrichment(uid, size_limit, raise_warning)
        cache.set(f'{uid}/analysis_enrichment', enrichment)

    info = dict(enrichment['info'])

    fieldnames = ['tf1', 'query1', 'analysis1', 'count1',
                  'tf2', 'query2', 'analysis2', 'count2',
                  'intersection_count', 'less', 'less_adj', 'greater', 'greater_adj', 'genes']

    writer = csv.DictWriter(buffer,
                            dialect='unix',
                            quoting=csv.QUOTE_MINIMAL,
                            fieldnames=fieldnames)

    writer.writeheader()

    for p, d, g in zip(enrichment['columns'],
                       map(itemgetter('less', 'less_adj', 'greater', 'greater_adj'), enrichment['data']),
                       map(lambda x: (len(x['genes']), ','.join(x['genes'])), enrichment['data'])):
        writer.writerow(dict(zip(fieldnames, p[0] + (info[p[0]]['Count'],) + p[1] + (info[p[1]]['Count'],) + d + g)))

    return buffer
