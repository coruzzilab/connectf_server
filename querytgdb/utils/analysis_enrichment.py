import csv
import warnings
from collections import OrderedDict
from functools import reduce
from io import StringIO
from itertools import combinations
from operator import itemgetter, methodcaller, or_
from typing import Dict, List, Optional, Tuple, Union
from uuid import UUID

import pandas as pd
from django.core.cache import cache
from django.http import HttpResponse
from scipy.special import comb
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

from querytgdb.utils import async_loader
from ..utils import clear_data, get_metadata


class AnalysisEnrichmentError(ValueError):
    pass


class AnalysisEnrichmentWarning(AnalysisEnrichmentError, UserWarning):
    pass


def split_col_name(col_name: Tuple[Tuple[str, str, str], int]) -> Tuple[str, str, str, int]:
    return col_name[0] + (col_name[1],)


def analysis_enrichment(uid: Union[UUID, str], size_limit: int = 100, raise_warning: bool = False) -> Dict:
    try:
        cached_data = cache.get_many([f'{uid}/tabular_output', f'{uid}/analysis_ids'])
        df, ids = itemgetter(
            f'{uid}/tabular_output',
            f'{uid}/analysis_ids'
        )(cached_data)
    except KeyError as e:
        raise AnalysisEnrichmentError("Please make a new query") from e

    df = df.pipe(clear_data)

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

    background_genes = cache.get(f'{uid}/background_genes')

    if background_genes is not None:
        background = background_genes.size
    else:
        background = async_loader['annotations'].shape[0]

    metadata = get_metadata(df.columns.get_level_values(1))

    for col_name in df.columns:
        d = OrderedDict(Count=df[col_name].count())
        d['filter'] = col_name[0][1]
        try:
            rename = ids[col_name]
            d['label'] = rename['name'] if rename['version'] else ''
        except KeyError as e:
            raise AnalysisEnrichmentError("Analysis Id data does not exist") from e

        d.update(metadata.loc[col_name[1], :].to_dict())

        info.append((split_col_name(col_name), d))

    for (name1, col1), (name2, col2) in combinations(((name, col.index[col.notna()])
                                                      for name, col in df.items()), 2):
        columns.append((
            split_col_name(name1),
            split_col_name(name2)
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
                .apply(lambda x: multipletests(x, method='bonferroni')[1])
                .rename(columns=lambda x: x + '_adj'))

    for d, adj_p in zip(data, p_values.itertuples(index=False)):
        d.update(adj_p._asdict())

    return {
        'columns': columns,
        'data': data,
        'info': info
    }


get_name_fields = itemgetter(0, 1, 3)
FIELD_NAMES = ['tf1', 'query1', 'analysis1', 'count1',
               'tf2', 'query2', 'analysis2', 'count2',
               'less', 'less_adj', 'greater', 'greater_adj', 'intersection_count', 'genes']


def analysis_enrichment_csv(uid: Union[str, UUID],
                            fields: Optional[List[str]] = None,
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

    if fields is not None:
        known_fields = reduce(or_, map(methodcaller('keys'), info.values()))
        fields = [f for f in fields if f in known_fields]
        fieldnames = FIELD_NAMES + [f + n for n in ('_1', '_2') for f in fields]
    else:
        fieldnames = FIELD_NAMES

    writer = csv.DictWriter(buffer,
                            dialect='unix',
                            quoting=csv.QUOTE_MINIMAL,
                            fieldnames=fieldnames)

    writer.writeheader()

    for p, d, g in zip(enrichment['columns'],
                       map(itemgetter('less', 'less_adj', 'greater', 'greater_adj'), enrichment['data']),
                       map(lambda x: (len(x['genes']), ','.join(x['genes'])), enrichment['data'])):
        row = get_name_fields(p[0]) + (info[p[0]]['Count'],) + get_name_fields(p[1]) + (info[p[1]]['Count'],) + d + g

        if fields is not None:
            row += tuple(info[p[n]][f] for n in (0, 1) for f in fields)

        writer.writerow(dict(zip(fieldnames, row)))

    return buffer
