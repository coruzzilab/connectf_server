import re
from functools import lru_cache
from itertools import groupby, islice, takewhile, zip_longest
from operator import itemgetter
from typing import Any, Dict, Iterable, List, Tuple, Union

import numpy as np
import pandas as pd

from ..models import Analysis


def is_numeric_column(cols: np.ndarray) -> np.ndarray:
    return np.isin(cols, ['User List Count', 'TF Count', 'Pvalue', 'Log2FC'])


def is_p_value(col: str) -> bool:
    return col == 'Pvalue'


def is_edge(col: Union[str, Tuple[str]]) -> bool:
    return col in {'EDGE', 'Log2FC'}


def is_data_column(col: str) -> bool:
    return col in {'Full Name', 'Gene Family', 'Type', 'Name', 'User List', 'User List Count', 'TF Count'}


def get_data_columns(df: pd.DataFrame) -> List[str]:
    return list(takewhile(is_data_column, df.columns))


def get_merge_cells(columns: List[Iterable]) -> List[Dict[str, Any]]:
    """
    Get merge cells for Handsontable. Highly customized don't copy paste.
    :param columns:
    :return:
    """
    columns = list(zip(*columns))

    merged_cells = [{'row': 0, 'col': i, 'colspan': 1, 'rowspan': 6} for i in range(8)]

    for i in range(1, 6):
        index = 8
        for label, group in groupby(columns[8:], key=itemgetter(slice(1, i + 1))):
            size = sum(1 for _ in group)

            if size > 1:
                merged_cells.append({'row': i, 'col': index, 'colspan': size, 'rowspan': 1})
                if i == 1:
                    merged_cells.append({'row': 0, 'col': index, 'colspan': size, 'rowspan': 1})

            index += size

    return merged_cells


def induce_repress_count(s: pd.Series) -> pd.Series:
    return pd.Series(((s > 0).sum(), (s < 0).sum()),
                     index=['induced', 'repressed'])


def get_col_type(s: Union[Tuple, str]) -> str:
    if isinstance(s, tuple):
        return s[2]
    else:
        return s


@lru_cache()
def get_tech(a: Analysis) -> str:
    return 'Data: {0.tech} (ID: {0.pk})'.format(a)


@lru_cache()
def get_analysis_method(a: Analysis) -> str:
    return 'Analysis: ' + a.analysis_method


@lru_cache()
def get_edge(a: Analysis) -> str:
    d = a.meta_dict

    return 'Edges: {0[EDGE_TYPE]} {{}}'.format(d)


def format_data(df: pd.DataFrame, stats: Dict) -> Tuple[List, List, List]:
    data_cols = get_data_columns(df)

    df = df.reset_index()
    df = df.reindex([*data_cols, df.columns[0], *df.columns[len(data_cols) + 1:]], axis=1)
    df = df.rename(columns={'TARGET': 'Gene ID'})

    col_types = np.fromiter(map(get_col_type, df.columns), 'U20', df.shape[1])

    num_cols = is_numeric_column(col_types)
    p_values = col_types == 'Pvalue'

    edge_cols = np.isin(col_types, ['EDGE', 'Log2FC'])
    edge_counts = df.loc[:, edge_cols].rename(columns=itemgetter(slice(None, 2))).count(axis=0)

    total_edge_counts = stats['total']

    fc_cols = col_types == 'Log2FC'
    induce_repress = df.loc[:, fc_cols].apply(induce_repress_count)

    empty_cols = df.isnull().all(axis=0)

    # for JSON response, can't have NaN or Inf
    df.loc[:, num_cols] = df.loc[:, num_cols].mask(np.isinf(df.loc[:, num_cols]), None)
    df = df.where(pd.notna(df), None)

    data_col_len = len(data_cols) + 1

    columns = list(map(list, zip_longest(*((col,) for col in df.columns[:data_col_len]),
                                         *df.columns[data_col_len:])))

    none_cols = [None] * len(columns[0])

    columns[-1: -1] = [none_cols.copy(),
                       none_cols.copy(),
                       none_cols.copy()]  # avoid references

    analyses = Analysis.objects.filter(
        pk__in=(c for c in columns[1] if c is not None)
    ).prefetch_related('analysisdata_set', 'analysisdata_set__key', 'tf')

    prev = None
    edge = None
    ind_rep = None
    for i, col in enumerate(islice(zip(*columns), data_col_len, None)):
        j = i + data_col_len
        name, _, uuid_ = col[0].rpartition(' ')
        name = name or uuid_

        try:
            analysis = analyses.get(pk=col[1])

            tf = analysis.tf
            if tf.name:
                name = re.sub(r'^' + re.escape(tf.gene_id), "{0.gene_id}<br/>({0.name})<br/>".format(tf), name, re.I)
            columns[0][j] = name
            columns[1][j] = get_tech(analysis)
            columns[2][j] = get_analysis_method(analysis)

            if col != prev:
                prev = col

                edge_count = edge_counts[col[:2]]
                total_edge_count = total_edge_counts[col[:2]]

                edge = get_edge(analysis).format(edge_count)
                if edge_count != total_edge_count:
                    edge += " ({})".format(total_edge_count)

                try:
                    ind_rep = "Induced-{0[induced]:} Repressed-{0[repressed]:}".format(
                        induce_repress[(*col[:2], 'Log2FC')])
                except KeyError:
                    ind_rep = None
            columns[3][j] = edge
            columns[4][j] = ind_rep
        except Analysis.DoesNotExist:
            pass

        if col[-1] == 'DAP':
            columns[2][j] = columns[1][i] = None

    merged_cells = get_merge_cells(columns)

    # Column formatting for Handsontable
    column_formats = []
    for i, (num, p, fc, empty, col) in enumerate(zip(num_cols, p_values, fc_cols, empty_cols, zip(*columns))):
        opt = {}
        if num:
            if p:
                opt['type'] = 'p_value'
            else:
                opt['type'] = 'numeric'

            if fc:
                opt['renderer'] = 'renderFc'
            elif i >= data_col_len and not p:
                opt['renderer'] = 'renderExp'

            if i >= data_col_len:
                opt['validator'] = 'exponential'
        else:
            opt['type'] = 'text'

            if col[-1] == 'EDGE':
                opt['className'] = 'htCenter'

        if empty:
            opt['width'] = 1

        column_formats.append(opt)

    return column_formats, merged_cells, [*columns, *df.itertuples(index=False, name=None)]
