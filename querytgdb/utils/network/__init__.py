import math
import os
from io import BytesIO
from operator import methodcaller
from typing import Any, BinaryIO, Dict, Generator, Iterable, List, Optional, Sized, Tuple
from uuid import uuid4

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import auc

from ...models import Analysis, Annotation
from ...utils import cache_result, data_to_edges, get_size, read_cached_result
from ...utils.network.utils import COLOR, COLOR_SHAPE
from ...utils.parser import ANNOTATIONS

GENE_TYPE = ANNOTATIONS[['Name', 'Type']]
SIZE = 20
GAP = 10
TF_GAP = 50
G_GAP = 40


def make_nodes(df: pd.DataFrame, pos, show_label: bool = False) -> Generator[Dict[str, Any], None, None]:
    """
    Make nodes for cytoscape network

    :param df:
    :param pos:
    :param show_label:
    :return:
    """
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


def group_edge_len(n: int, size: int = SIZE, gap: int = GAP) -> int:
    return n * size + (n - 1) * gap


def get_network_json(df: pd.DataFrame, edges: Optional[List[str]] = None) -> List[Dict[str, Any]]:
    """
    Get cytoscape network from queried data.

    :param df:
    :param edges:
    :return:
    """
    network_table = df.loc[:, (slice(None), slice(None), ['EDGE', 'Log2FC'])]

    analyses = Analysis.objects.filter(
        pk__in=df.columns.get_level_values(1)
    ).prefetch_related('analysisdata_set', 'analysisdata_set__key', 'tf')

    network_table = data_to_edges(network_table, analyses)

    def get_tf(idx):
        return analyses.get(pk=idx).tf.gene_id

    network_table = network_table.rename(columns=get_tf, level=1)
    network_table.columns = network_table.columns.droplevel(0)

    # add additional edges here

    edge_nodes = network_table.reindex(index=network_table.index.difference(network_table.columns))

    uniq_tfs = network_table.columns.unique()

    try:
        tf_nodes = network_table.reindex(index=uniq_tfs)
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

    data.extend({
                    'group': 'edges',
                    'data': {
                        'id': uuid4(),
                        'source': s,
                        'target': t,
                        'name': e,
                        'color': COLOR[e]
                    }
                } for t, s, e in tf_nodes.stack().reset_index().drop_duplicates().itertuples(name=None, index=False))

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

        data.extend({
                        'group': 'edges',
                        'data': {
                            'id': uuid4(),
                            'source': s,
                            'target': t,
                            'name': e,
                            'color': COLOR[e],
                            'weight': w
                        }
                    } for t, s, e, w in e_group.reset_index().itertuples(name=None, index=False))

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


def get_network_stats(df: pd.DataFrame) -> Dict[str, Any]:
    df = df.loc[:, (slice(None), slice(None), ['EDGE', 'Log2FC'])]

    tfs = Annotation.objects.filter(
        analysis__in=df.columns.get_level_values(1)
    ).distinct().count()

    return {
        'num_edges': df.count().sum(),
        'num_targets': df.shape[0],
        'num_tfs': tfs
    }


def fix_tied(dup: np.ndarray, s: pd.Series) -> pd.Series:
    """
    Replace all duplicates (tied ranks) with last value
    :param dup:
    :param s:
    :return:
    """
    return s.mask(dup).fillna(method='bfill')


AucData = Tuple[float, Sized, Sized]
row_labels = [
    'AUPR',
    'AUPR Random',
    'p-value',
    'Precision Threshold',
    'Edge Score Threshold',
    'Number of edges',
    'Number of TFs',
    'Number of Targets'
]


def get_auc(network: Tuple[str, pd.DataFrame], df: pd.DataFrame, precision_cutoff: Optional[float] = None,
            cache_path: Optional[str] = None) -> BinaryIO:
    """
    Get AUPR of uploaded predicted network

    :param network:
    :param df:
    :param precision_cutoff:
    :param cache_path:
    :return:
    """
    fig: plt.Figure
    gs: plt.GridSpec
    recall: pd.Series
    precision: pd.Series
    g: pd.DataFrame
    cell_text: List[List[str]]

    buff = BytesIO()  # place to store the actual graph
    name, data = network
    data = data.sort_values('rank')

    try:
        cache_file = os.path.join(cache_path, 'figure.pickle.gz')

        fig, gs, recall, precision, g, cell_text = read_cached_result(cache_file)
        plt.figure(fig.number)

    except (TypeError, FileNotFoundError) as e:
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

        cell_text = [[''] for i in range(8)]

        fig = plt.figure(figsize=(9.6, 4.8))
        gs = plt.GridSpec(1, 2, width_ratios=[1, 2])

        plt.subplot(gs[1])

        plt.xlabel("recall")
        plt.xlim(-0.05, 1.05)
        plt.ylabel("precision")
        plt.ylim(-0.05, 1.05)

        d = df.loc[df['TARGET'].isin(data['target']) & df['ANALYSIS'].isin(data['source']), :]
        g = data.loc[data['source'].isin(d['ANALYSIS']), :]  # filtering predictions
        g = g.merge(d, how='left', left_on=['source', 'target'], right_on=['ANALYSIS', 'TARGET']).sort_values('rank')
        g[0] = g[0].fillna(0).astype(int)

        dup = g['rank'].duplicated(keep='last').values  # use ndarray to disregard indices

        r = np.arange(1, g.shape[0] + 1)
        c = g[0].cumsum()
        recall = fix_tied(dup, c / g[0].sum())
        precision = fix_tied(dup, c / r)

        pred_auc = auc(recall, precision)

        rand_aucs = []

        h = g[0]

        min_auc: AucData = (1, [], [])
        max_auc: AucData = (0, [], [])

        for i in range(1000):
            h = h.sample(frac=1)  # shuffle
            c = h.cumsum()
            curr_prec = fix_tied(dup, c / r)
            curr_recall = fix_tied(dup, c / h.sum())
            curr_auc = auc(curr_recall, curr_prec)
            rand_aucs.append(curr_auc)

            if curr_auc <= min_auc[0]:
                min_auc = (curr_auc, curr_recall, curr_prec)

            if curr_auc >= max_auc[0]:
                max_auc = (curr_auc, curr_recall, curr_prec)

        p_value = (np.array(rand_aucs) >= pred_auc).mean()

        plt.plot(recall, precision, label=f'{name} auc: {pred_auc:.4f}', zorder=3, color='C0')
        plt.plot(*max_auc[1:], label=f'{name} max random auc: {max_auc[0]:.4f}',
                 linestyle='--', zorder=1, color='darkgrey')
        plt.plot(*min_auc[1:], label=f'{name} min random auc: {min_auc[0]:.4f}',
                 linestyle=':', zorder=2, color='lightgrey')

        # setup the table style
        ax = plt.subplot(gs[0])
        ax.patch.set_visible(False)
        ax.axis('off')

        # Coordinate with precision_cutoff
        cell_text[:3] = [
            [format(pred_auc, '.4f')],
            [format(np.mean(rand_aucs), '.4f')],
            [format(p_value, '.3f') if p_value else '<0.001']
        ]

        try:
            cache_result((fig, gs, recall, precision, g, cell_text), cache_file)
        except NameError:
            pass

    if precision_cutoff is not None:
        cell_text[3][0] = str(precision_cutoff)

        plt.subplot(gs[1])
        plt.axhline(y=precision_cutoff, color='red', label=f"precision cutoff: {precision_cutoff}", linestyle='--')
        # plt.annotate(precision_cutoff, (1, precision_cutoff), color='red')

        cutoff_loc: pd.Series = (precision >= precision_cutoff)

        if cutoff_loc.any():
            max_idx = cutoff_loc.iloc[::-1].idxmax()
            max_loc = g.index.get_loc(max_idx)

            xy = recall.iloc[max_loc], precision.iloc[max_loc]

            plt.plot(*xy, 'ro', fillstyle='none')

            rank = g.at[max_idx, "rank"]
            rank_loc = data["rank"] <= rank

            cell_text[5:] = [
                ["{:,}/{:,}".format(rank_loc.sum(), data.shape[0])],
                ["{:,}/{:,}".format(data.loc[rank_loc, "source"].nunique(), data["source"].nunique())],
                ["{:,}/{:,}".format(data.loc[rank_loc, "target"].nunique(), data["target"].nunique())]
            ]

            s = f'precision: {xy[1]:.04}\nrecall: {xy[0]:0.4}'

            try:
                score = format(g.at[max_idx, "score"], '.4f')
                cell_text[4][0] = score
                s += f'\nscore: {score}'
            except KeyError:
                pass

            plt.annotate(s, xy, xytext=(3, 3), textcoords='offset pixels', color='red')

        else:
            cell_text[5:] = [
                ["0/{:,}".format(data.shape[0])],
                ["0/{:,}".format(data["source"].nunique())],
                ["0/{:,}".format(data["target"].nunique())]
            ]
    else:
        cell_text[5:] = [
            ["{:,}".format(data.shape[0])],
            ["{:,}".format(data["source"].nunique())],
            ["{:,}".format(data["target"].nunique())]
        ]

    plt.subplot(gs[0])
    plt.table(
        cellText=cell_text,
        rowLabels=row_labels,
        loc='center'
    )

    plt.subplot(gs[1])
    plt.legend(loc=0)
    plt.savefig(buff, bbox_inches='tight')

    plt.close(fig)
    buff.seek(0)

    return buff
