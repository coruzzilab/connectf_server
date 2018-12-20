import math
import re
from collections import UserDict
from io import BytesIO
from operator import methodcaller
from typing import Any, BinaryIO, Dict, Generator, Iterable, List, Union
from uuid import uuid4

import hsluv
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import auc

from querytgdb.utils import get_size
from ..models import Analysis
from ..utils import data_to_edges
from ..utils.parser import ANNOTATIONS
from ..utils.file import Graphs


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

COLOR = {
    'INDUCED': '#4daf4a',
    'REPRESSED': '#e41a1c',
    'BOUND': '#377eb8'
}


def simplify_edge(edge: str) -> str:
    if edge.endswith("INDUCED"):
        return "INDUCED"
    elif edge.endswith("REPRESSED"):
        return "REPRESSED"
    return "BOUND"


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


def get_network_json(df: pd.DataFrame, graphs: Graphs = None) -> List[Dict[str, Any]]:
    network_table = df.loc[:, (slice(None), slice(None), ['EDGE', 'Log2FC'])]

    network_table = data_to_edges(network_table)

    network_table.columns = network_table.columns.droplevel(1)
    network_table = network_table.rename(columns=lambda x: re.split(r'\s', x, 1)[0], level=0)

    edge_nodes = network_table.loc[network_table.index.difference(network_table.columns.get_level_values(0)), :]

    uniq_tfs = network_table.columns.get_level_values(0).unique()

    try:
        tf_nodes = network_table.loc[uniq_tfs, :]
    except KeyError:
        tf_nodes = pd.DataFrame(index=uniq_tfs, columns=network_table.columns)
    tf_nodes = tf_nodes.reindex(index=tf_nodes.count(axis=1).sort_values(ascending=False).index)

    # having duplicate names is a no-no so clear all names
    tf_nodes.columns.name = None
    tf_nodes.index.name = None

    s_tfs = math.ceil(math.sqrt(tf_nodes.shape[0]))
    e_tfs = group_edge_len(s_tfs, gap=TF_GAP)

    tf_grid = np.array(np.meshgrid(np.arange(s_tfs), np.arange(s_tfs), indexing='ij')).reshape(
        (2, -1), order='F').T * (SIZE + TF_GAP) + SIZE / 2

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
        'color': COLOR[simplify_edge(e)]
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
            'color': COLOR[simplify_edge(e)],
            'weight': w
        }} for t, s, e, w in e_group.reset_index().itertuples(name=None, index=False))

    return data


@get_size
def get_points(recall: Iterable[float], precision: Iterable[float]) -> List[Dict[str, float]]:
    """
    Format data points for the json output
    :param recall:
    :param precision:
    :return:
    """
    return [{'x': r, 'y': p} for r, p in zip(recall, precision)]


def get_auc(graphs: Graphs, df: pd.DataFrame) -> Union[BinaryIO, Dict]:
    df = df.loc[:, (slice(None), slice(None), ['EDGE', 'Log2FC'])]
    df.columns = df.columns.droplevel(2)
    df = df.notna()

    analyses = Analysis.objects.filter(pk__in=df.columns.get_level_values(1)).prefetch_related('tf')

    def get_tf(pk):
        return analyses.get(pk=pk).tf.gene_id

    df = df.groupby(get_tf, axis=1, level=1).apply(methodcaller('any', axis=1))
    df.columns = df.columns.str.upper()
    df.index = df.index.str.upper()
    df[~df] = np.nan  # query data matrix
    df = df.stack().reset_index()
    df.columns = ['TARGET', 'ANALYSIS', 0]

    # result = {}

    buff = BytesIO()

    plt.xlabel("recall")
    plt.xlim(-0.05, 1.05)
    plt.ylabel("precision")
    plt.ylim(-0.05, 1.05)

    for name, graph in graphs.items():
        g = pd.DataFrame(iter(graph.edges.data('rank')))
        d = df.loc[df['TARGET'].isin(g[1]) & df['ANALYSIS'].isin(g[0]), :]
        g = g.loc[g[0].isin(d['ANALYSIS']), :]  # filtering predictions
        g = g.merge(d, how='left', left_on=[0, 1], right_on=['ANALYSIS', 'TARGET']).sort_values(2)
        g['0_y'] = g['0_y'].fillna(0).astype(int)

        r = (np.arange(g.shape[0]) + 1)
        recall = g['0_y'].cumsum() / g['0_y'].sum()
        precision = g['0_y'].cumsum() / r

        # max_loc = g.index.get_loc((precision > 0.7).iloc[::-1].idxmax())

        # g.iloc[:max_loc + 1, :]

        pred_auc = auc(recall, precision)

        plt.plot(recall, precision, label=f'{name}  auc: {pred_auc:.4f}', zorder=3)

        rand_aucs = []

        h = g['0_y']

        min_auc = (1, [], [])
        max_auc = (0, [], [])

        for i in range(1000):
            h = h.sample(frac=1)  # shuffle
            c = h.cumsum()
            curr_prec = c / r
            curr_recall = c / h.sum()
            curr_auc = auc(curr_recall, curr_prec)
            rand_aucs.append(curr_auc)

            if curr_auc <= min_auc[0]:
                min_auc = (curr_auc, curr_recall, curr_prec)

            if curr_auc >= max_auc[0]:
                max_auc = (curr_auc, curr_recall, curr_prec)

        # result[name] = {
        #     'auc': pred_auc,
        #     'data': get_points(recall, precision),
        #     'min_auc': (min_auc[0], get_points(*min_auc[1:])),
        #     'max_auc': (max_auc[0], get_points(*max_auc[1:])),
        #     'p_value': (np.array(rand_aucs) > pred_auc).mean()
        # }

        plt.plot(max_auc[1], max_auc[2], label=f'{name} max random auc: {max_auc[0]:.4f}', linestyle='--', zorder=1)
        plt.plot(min_auc[1], min_auc[2], label=f'{name} min random auc: {min_auc[0]:.4f}', linestyle=':', zorder=2)

    plt.legend()
    plt.savefig(buff)
    plt.close()
    buff.seek(0)

    return buff
