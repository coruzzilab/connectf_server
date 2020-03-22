import logging
import re
from itertools import groupby, islice, zip_longest
from operator import itemgetter
from typing import Any, Dict, Iterable, List, Tuple, Union

import numpy as np
import pandas as pd

from ..utils.parser import Ids

logger = logging.getLogger(__name__)


def is_numeric_column(cols: np.ndarray) -> np.ndarray:
    return np.isin(cols, ['User List Count', 'TF Count', 'Pvalue', 'Log2FC'])


def is_p_value(col: str) -> bool:
    return col == 'Pvalue'


def is_edge(col: Union[str, Tuple[str]]) -> bool:
    return col in {'EDGE', 'Log2FC'}


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


def get_col_type(s: Union[Tuple, str]) -> str:
    if isinstance(s, tuple):
        return s[2]
    else:
        return s


def get_name(name, a: pd.Series, rename: Dict[str, Any]) -> str:
    name = re.sub(r'^' + re.escape(a['gene_id']), a['gene_id'], name, 1, re.I)
    if a['gene_name'] and rename['version']:
        name += " ({0[gene_name]})\n{1[name]}\n".format(a, rename)
    elif a['gene_name']:
        name += "\n({0[gene_name]})\n".format(a)
    elif rename['version']:
        name += "\n{[name]}\n".format(rename)

    return name


def get_tech(a: pd.Series) -> str:
    return 'Data: {} (ID: {})'.format(
        a.get('EXPERIMENTER', default='') + a.get('DATE', default='').replace('-', '') + '_' + a.get('TECHNOLOGY',
                                                                                                     default=''),
        a.get('analysis_id', default='')
    )


def get_analysis_method(a: pd.Series) -> str:
    return 'Analysis: ' + a.get('ANALYSIS_METHOD', default='') + '_' + re.sub(r'\s+', '',
                                                                              a.get('ANALYSIS_CUTOFF', default=''))


def get_edge(a: pd.Series) -> str:
    return 'Edges: {0} {{}}'.format(a.get('EDGE_TYPE', default=''))


def format_data(df: pd.DataFrame, stats: Dict, metadata: pd.DataFrame, ids: Ids) -> Tuple[List, List, List]:
    df = df.reset_index()
    df.insert(7, 'Gene ID', df.pop('TARGET'))

    col_types = np.fromiter(map(get_col_type, df.columns), 'U20', df.shape[1])

    num_cols = is_numeric_column(col_types)
    num_df = df.loc[:, num_cols]
    p_values = col_types == 'Pvalue'

    edge_counts = stats['edge_counts']
    edge_counts.index = edge_counts.index.tolist()

    total_edge_counts = stats['total']
    total_edge_counts.index = total_edge_counts.index.tolist()

    fc_cols = col_types == 'Log2FC'
    induce_repress = stats['induce_repress_count']
    induce_repress.index = induce_repress.index.tolist()

    empty_cols = df.isnull().all(axis=0)

    # for JSON response, can't have NaN or Inf
    df.loc[:, num_cols] = num_df.mask(np.isinf(num_df), None)
    df = df.where(pd.notna(df), None)

    data_col_len = 8

    columns = list(map(list, zip_longest(*((col,) for col in df.columns[:data_col_len]),
                                         *df.columns[data_col_len:])))

    none_cols = [None] * len(columns[0])

    columns[-1: -1] = [none_cols.copy(),
                       none_cols.copy(),
                       none_cols.copy()]  # avoid references

    prev = (None,) * 5
    name = None
    tech = None
    method = None
    edge = None
    ind_rep = None
    for i, col in enumerate(islice(zip(*columns), data_col_len, None), data_col_len):
        try:
            if col[0:2] != prev[0:2]:
                name, _, uuid_ = col[0]
                analysis = metadata.loc[:, col[1]]

                name = get_name(name, analysis, ids[col[0:2]])

                tech = get_tech(analysis)
                method = get_analysis_method(analysis)

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

            columns[0][i] = name
            columns[1][i] = tech
            columns[2][i] = method
            columns[3][i] = edge
            columns[4][i] = ind_rep
        except KeyError:
            pass

        if col[-1] == 'DAP':
            columns[2][i] = columns[1][i] = None
        prev = col

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
                opt['renderer'] = 'renderEdge'

        if empty:
            opt['width'] = 1

        column_formats.append(opt)
    return column_formats, merged_cells, columns + df.values.tolist()
