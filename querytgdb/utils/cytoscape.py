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
TF_GAP = 50
G_GAP = 40


def make_nodes(df: pd.DataFrame, pos, show_label: bool = False) -> Generator[Dict[str, Any], None, None]:
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
    network_table = network_table.rename(columns=lambda x: re.split(r'\s', x, 1)[0], level=0)

    edge_nodes = network_table.loc[network_table.index.difference(network_table.columns.get_level_values(0)), :]

    uniq_tfs = network_table.columns.get_level_values(0).unique()
    try:
        tf_nodes = network_table.loc[uniq_tfs, :]
    except KeyError:
        tf_nodes = pd.DataFrame(index=uniq_tfs, columns=network_table.columns)
    tf_nodes = tf_nodes.loc[tf_nodes.count(axis=1).sort_values(ascending=False).index, :]

    s_tfs = math.ceil(math.sqrt(tf_nodes.shape[0]))
    e_tfs = group_edge_len(s_tfs, gap=TF_GAP)

    tf_grid = np.array(np.meshgrid(np.arange(s_tfs), np.arange(s_tfs), indexing='ij')).reshape(
        (2, -1), order='F').T * (SIZE + TF_GAP) + SIZE / 2

    edge_colors = get_edge_colors(network_table)

    edge_node_stack = edge_nodes.stack()
    # edge building and positioning

    edge_counts = (edge_node_stack
                   .groupby(level=[0, 1])
                   .value_counts())
    # scale edge weights here
    edge_counts = edge_counts - edge_counts.min() + 1

    edge_group = (edge_counts
                  .reset_index(level=[1, 2])
                  .groupby(edge_node_stack
                           .reset_index(level=1)
                           .groupby(level=0)
                           .apply(lambda x: x.iloc[:, 0].nunique())))

    max_group_num = edge_group.apply(lambda x: x.index.nunique()).max()

    num_targets = len(edge_group)
    s_target = math.ceil(math.sqrt(num_targets))

    num_group = math.ceil(math.sqrt(max_group_num))
    group_bbox = group_edge_len(num_group)  # square edge length of largest number of tfs

    groups_edge_len = group_edge_len(s_target, group_bbox, G_GAP)

    tf_grid += ((groups_edge_len - e_tfs) / 2, 0)
    data = list(make_nodes(GENE_TYPE.loc[tf_nodes.index], tf_grid, True))

    data.extend({'group': 'edges', 'data': {
        'id': uuid4(),
        'source': s,
        'target': t,
        'name': e,
        'color': edge_colors[e]
    }} for t, s, e in tf_nodes.stack().reset_index().drop_duplicates().itertuples(name=None, index=False))

    group_grid = np.array(np.meshgrid(np.arange(s_target), np.arange(s_target), indexing='ij')).reshape(
        (2, -1), order='F').T * (group_bbox + G_GAP) + group_bbox / 2 + (0, (groups_edge_len + e_tfs) / 2)

    for (num, e_group), g_ij in zip(edge_group, group_grid):
        e_group = e_group.sort_values(by=e_group.columns.tolist())
        s_group = math.ceil(math.sqrt(e_group.index.nunique()))

        e_grid = np.array(np.meshgrid(np.arange(s_group), np.arange(s_group), indexing='ij'), dtype=np.float64).reshape(
            (2, -1), order='C').T * (SIZE + GAP)
        e_grid += g_ij - np.mean(e_grid, axis=0)  # add offset

        data.extend(
            make_nodes(
                GENE_TYPE.loc[e_group.index.unique(), :].sort_values('Type', kind='mergesort'),
                e_grid))

        data.extend({'group': 'edges', 'data': {
            'id': uuid4(),
            'source': s,
            'target': t,
            'name': e,
            'color': edge_colors[e],
            'weight': w
        }} for t, s, e, w in e_group.reset_index().itertuples(name=None, index=False))

    return data
