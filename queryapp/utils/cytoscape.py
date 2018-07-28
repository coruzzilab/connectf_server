import math
import re
from collections import UserDict
from typing import Any, Dict, Generator, List
from uuid import uuid4

import hsluv
import numpy as np
import pandas as pd

from ..utils.parser import ANNOTATIONS


class Shape(UserDict):
    def __init__(self):
        super().__init__({
            'TXNFACTOR': {
                'color': '#00FF00',
                'shape': 'triangle'
            },
            'PROTEIN_CODING': {
                'color': '#AED6F1',
                'shape': 'roundrectangle'
            },
            'METABOLIC': {
                'color': '#D0ECE7',
                'shape': 'roundrectangle'
            },
            'MOLECULE': {
                'color': '#FF9900',
                'shape': 'roundrectangle'
            }
        })

    def __getitem__(self, item):
        try:
            return super().__getitem__(item)
        except KeyError:
            return {
                'color': '#FF9900',
                'shape': 'roundrectangle'
            }


GENE_TYPE = ANNOTATIONS[['Name', 'Type']]
COLOR_SHAPE = Shape()
SIZE = 20
GAP = 10
G_GAP = 40


def make_nodes(df: pd.DataFrame, pos, show_label: str = 'hide') -> Generator[Dict[str, Any], None, None]:
    for (idx, name, t), (y, x) in zip(df.itertuples(name=None), pos):
        yield {
            'group': 'nodes',
            'data': {
                'id': idx,
                'name': name,
                'type': t,
                'size': SIZE,
                'showLabel': show_label,
                **COLOR_SHAPE[t]
            },
            'position': {
                'x': x,
                'y': y
            }
        }


def get_edge_colors(df: pd.DataFrame, s: float = 100, l: float = 50) -> Dict[str, str]:
    edges = df.values.ravel('K')
    edges = pd.unique(edges[pd.notna(edges)])

    return {e: hsluv.hsluv_to_hex((i, s, l)) for e, i in zip(edges, np.linspace(0, 360, len(edges), False))}


def group_edge_len(n: int, size: int = SIZE, gap: int = GAP) -> int:
    return n * size + (n - 1) * gap


def get_cytoscape_json(df: pd.DataFrame) -> List[Dict[str, Any]]:
    network_table = df.loc[:, (slice(None), slice(None), slice(None), 'EDGE')]
    network_table.columns = network_table.columns.droplevel(level=[1, 2, 3])
    network_table = network_table.rename(columns=lambda x: re.split(r'[\s_]', x, 1)[0], level=0)

    edge_nodes = network_table.loc[network_table.index.difference(network_table.columns.get_level_values(0)), :]

    uniq_tfs = network_table.columns.get_level_values(0).unique()
    try:
        tf_nodes = network_table.loc[uniq_tfs, :]
    except KeyError:
        tf_nodes = pd.DataFrame(index=uniq_tfs, columns=network_table.columns)
    tf_nodes = tf_nodes.loc[tf_nodes.count(axis=1).sort_values(ascending=False).index, :]

    s_tfs = math.ceil(math.sqrt(tf_nodes.shape[0]))
    e_tfs = group_edge_len(s_tfs)

    tf_grid = np.array(np.meshgrid(np.arange(s_tfs), np.arange(s_tfs), indexing='ij')).reshape(
        (2, -1), order='F').T * (SIZE + 50) + SIZE / 2

    edge_colors = get_edge_colors(network_table)

    # edge building and positioning
    edge_group = edge_nodes.groupby(edge_nodes.count(axis=1))

    max_group_num = edge_group.size().max()

    num_targets = len(edge_group)
    s_target = math.ceil(math.sqrt(num_targets))

    num_group = math.ceil(math.sqrt(max_group_num))
    group_bbox = group_edge_len(num_group)  # square edge length of largest number of tfs

    groups_edge_len = group_edge_len(s_target, group_bbox, G_GAP)

    tf_grid += ((groups_edge_len - e_tfs) / 2, 0)
    data = list(make_nodes(GENE_TYPE.loc[tf_nodes.index], tf_grid, 'show'))

    data.extend({'group': 'edges', 'data': {
        'id': uuid4(),
        'source': s,
        'target': t,
        'name': e,
        'color': edge_colors[e]
    }} for (t, s), e in tf_nodes.stack().iteritems())

    group_grid = np.array(np.meshgrid(np.arange(s_target), np.arange(s_target), indexing='ij')).reshape(
        (2, -1), order='F').T * (group_bbox + G_GAP) + group_bbox / 2 + (0, groups_edge_len / 4)

    for (num, e_group), g_ij in zip(edge_group, group_grid):
        stack_group = e_group.stack().sort_values()

        s_group = math.ceil(math.sqrt(e_group.shape[0]))

        e_grid = np.array(np.meshgrid(np.arange(s_group), np.arange(s_group), indexing='ij'), dtype=np.float64).reshape(
            (2, -1), order='C').T * (SIZE + GAP)
        e_grid += g_ij - np.mean(e_grid, axis=0)  # add offset

        data.extend(
            make_nodes(
                GENE_TYPE.loc[stack_group.index.unique(level=0), :].sort_values('Type', kind='mergesort'),
                e_grid))

        data.extend({'group': 'edges', 'data': {
            'id': uuid4(),
            'source': s,
            'target': t,
            'name': e,
            'color': edge_colors[e]
        }} for (t, s), e in stack_group.iteritems())

    return data
