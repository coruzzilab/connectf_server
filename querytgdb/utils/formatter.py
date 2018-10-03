from itertools import groupby, islice, takewhile, zip_longest
from operator import itemgetter
from typing import Dict, Iterable, List, Tuple, Union

import numpy as np
import pandas as pd


def is_numeric_column(col: str) -> bool:
    return col in {'User List Count', 'TF Count', 'Pvalue', 'Log2FC'}


def is_p_value(col: str) -> bool:
    return col == 'Pvalue'


def is_edge(col: Union[str, Tuple[str]]) -> bool:
    return col == 'EDGE'


def is_data_column(col: str) -> bool:
    return col in {'Full Name', 'Gene Family', 'Type', 'Name', 'User List', 'User List Count', 'TF Count'}


def get_data_columns(df: pd.DataFrame) -> List[str]:
    return list(takewhile(is_data_column, df.columns))


def get_merge_cells(columns: List[Iterable]) -> List:
    """
    Get merge cells for Handsontable. Highly customized don't copy paste.
    :param columns:
    :return:
    """
    columns = list(zip(*columns))

    merged_cells = [{'row': 0, 'col': i, 'colspan': 1, 'rowspan': 6} for i in range(8)]

    for i in range(6):
        index = 8
        for label, group in groupby(columns[8:], key=itemgetter(*range(0, i + 1))):
            size = sum(1 for _ in group)

            if size > 1:
                merged_cells.append({'row': i, 'col': index, 'colspan': size, 'rowspan': 1})

            index += size

    return merged_cells


def induce_repress_count(s: pd.Series):
    return pd.Series(((s > 0).sum(), (s < 0).sum()),
                     index=['induced', 'repressed'])


def get_col_type(s):
    if isinstance(s, tuple):
        return s[3]
    else:
        return s


def format_data(df: pd.DataFrame, stats: Union[Dict, None] = None):
    data_cols = get_data_columns(df)

    df = df.reset_index()
    df = df.reindex([*data_cols, df.columns[0], *df.columns[len(data_cols) + 1:]], axis=1)
    df = df.rename(columns={'TARGET': 'Gene ID'})

    col_types = list(map(get_col_type, df.columns))

    num_cols = list(map(is_numeric_column, col_types))
    p_values = list(map(is_p_value, col_types))

    edge_cols = list(map(is_edge, col_types))
    edge_counts = df.loc[:, edge_cols].count(axis=0)

    total_edge_counts = stats['total']

    fc_cols = list(map(lambda c: c == 'Log2FC', col_types))
    induce_repress = df.loc[:, fc_cols].apply(induce_repress_count)

    # for JSON response, can't have NaN or Inf
    df.loc[:, num_cols] = df.loc[:, num_cols].mask(np.isinf(df.loc[:, num_cols]), None)
    df = df.where(pd.notna(df), None)

    data_col_len = len(data_cols) + 1

    # Column formatting for Handsontable
    column_formats = []
    for i, (num, p, fc) in enumerate(zip(num_cols, p_values, fc_cols)):
        opt = {}
        if num:
            if p:
                opt.update({'type': 'p_value'})
            else:
                opt.update({'type': 'numeric'})

            if fc:
                opt.update({'renderer': 'renderFc'})
            elif i >= data_col_len and not p:
                opt.update({'renderer': 'renderExp'})

            if i >= data_col_len:
                opt.update({'validator': 'exponential'})
        else:
            opt.update({'type': 'text'})

        column_formats.append(opt)

    columns = list(map(list, zip_longest(*((col,) for col in df.columns[:data_col_len]),
                                         *df.columns[data_col_len:])))

    columns[-1: -1] = [[None] * len(columns[0]),
                       [None] * len(columns[0])]  # don't use "* 2" to make a copy instead of reference

    for i, col in enumerate(columns[-1]):
        if col == 'DAP':
            columns[2][i] = columns[1][i] = None

    prev = None
    edge = None
    ind_rep = None
    for i, col in enumerate(islice(zip(*columns), data_col_len, None)):
        name, _, uuid_ = col[0].rpartition(' ')
        columns[0][i + data_col_len] = name or uuid_
        if col != prev:
            prev = col

            edge_count = edge_counts[(*col[:3], 'EDGE')]
            total_edge_count = total_edge_counts[col[:3]]

            edge = "Edges: {}".format(edge_count)
            if edge_count != total_edge_count:
                edge += " ({})".format(total_edge_count)

            try:
                ind_rep = "Induced-{} Repressed-{}".format(*induce_repress[(*col[:3], 'Log2FC')])
            except KeyError:
                ind_rep = None
        columns[3][i + data_col_len] = edge
        columns[4][i + data_col_len] = ind_rep

    merged_cells = get_merge_cells(columns)

    return column_formats, merged_cells, [*columns, *df.itertuples(index=False, name=None)]
