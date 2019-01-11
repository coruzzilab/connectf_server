import sys
from itertools import count, product
from operator import itemgetter
from typing import Any, Dict, Optional
from uuid import uuid4

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as hierarchy
import seaborn as sns
from scipy.stats import fisher_exact

from ..models import Analysis
from ..utils import clear_data, column_string, split_name
from ..utils.parser import ANNOTATIONS


def scale_df(df: pd.DataFrame) -> pd.DataFrame:
    df = df.clip_lower(sys.float_info.min)
    return -np.log10(df)


def draw_heatmap(df: pd.DataFrame, **kwargs):
    row_num, col_num = df.shape
    opts = {}
    if row_num > 1:
        opts['row_linkage'] = hierarchy.linkage(df.values, method='average', optimal_ordering=True)
    if col_num > 1:
        opts['col_linkage'] = hierarchy.linkage(df.values.T, method='average', optimal_ordering=True)

    sns_heatmap = sns.clustermap(df,
                                 cmap="YlGnBu",
                                 cbar_kws={'label': 'Enrichment (-log10 p)'},
                                 xticklabels=1,
                                 yticklabels=1,
                                 row_cluster=row_num > 1,
                                 col_cluster=col_num > 1,
                                 **kwargs,
                                 **opts)
    plt.setp(sns_heatmap.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    # plt.setp(sns_heatmap.ax_heatmap.yaxis.get_label(), text="Analyses", rotation=270)
    plt.setp(sns_heatmap.ax_heatmap.xaxis.get_majorticklabels(), rotation=270)
    plt.setp(sns_heatmap.ax_heatmap.xaxis.get_label(), text="Target Genes")

    return sns_heatmap


def gene_list_enrichment(pickledir, background: Optional[int] = None, draw=True, legend=False, save_file=False,
                         upper=None, lower=None):
    """
    Draws gene list enrichment heatmap from cached query

    :param pickledir:
    :param background:
    :param draw:
    :param legend:
    :param save_file:
    :param upper:
    :param lower:
    :return:
    """
    # raising exception here if target genes are not uploaded by the user
    try:
        name_to_list, list_to_name = pd.read_pickle(pickledir + '/target_genes.pickle.gz')
    except FileNotFoundError as e:
        raise FileNotFoundError('No target genes uploaded') from e

    if background is None:
        background = ANNOTATIONS.shape[0]
        if not background:
            background = 28775

    query_result = pd.read_pickle(pickledir + '/tabular_output_unfiltered.pickle.gz')

    # default background:
    # A. thaliana columbia tair10 genome(28775 genes(does not include transposable elements and pseudogenes))

    query_result = clear_data(query_result)

    analyses = Analysis.objects.filter(
        pk__in=map(itemgetter(1), query_result.columns)
    ).distinct().prefetch_related('tf')

    # save targets for each TF into a dict
    targets = {}

    for (name, analysis_id), column in query_result.iteritems():
        analysis_targets = set(column.dropna().index)

        name, criterion = split_name(name)

        targets[(name, criterion, analysis_id, len(analysis_targets), uuid4())] = analysis_targets

    list_enrichment_pvals = pd.DataFrame(index=targets.keys(), columns=list_to_name.keys(), dtype=np.float64)

    if not legend:
        colnames = {}

        for (analysis_name, analysis_list), (name, user_list) in product(targets.items(), list_to_name.items()):
            intersect_len = len(analysis_list & user_list)
            colnames[name] = "{} ({})".format(name, len(user_list))

            odds, pvalue = fisher_exact([
                [intersect_len, len(user_list - analysis_list)],
                [len(analysis_list - user_list), background - len(user_list | analysis_list)]
            ], alternative='greater')
            list_enrichment_pvals.at[analysis_name, name] = pvalue

        list_enrichment_pvals.rename(columns=colnames, inplace=True)

    if draw or legend:
        orig_index = list(zip(list_enrichment_pvals.index, map(column_string, count(1))))

        indices = []

        for (name, criterion, analysis_id, l, uid), col_name in orig_index:
            tf = analyses.get(pk=analysis_id).tf
            indices.append('{} — {}{} ({})'.format(col_name, tf.gene_id, ' ' + tf.name if tf.name else '', l))

        list_enrichment_pvals.index = indices

        if legend:
            result = []

            for idx, col_label in orig_index:
                name, criterion, analysis_id, l, uid = idx

                info = {'name': name, 'filter': criterion}
                analysis = analyses.get(pk=analysis_id)

                gene_id = analysis.tf.gene_id

                try:
                    gene_name = ANNOTATIONS.at[gene_id, 'Name']
                except KeyError:
                    gene_name = ''

                info.update(analysis.analysisdata_set.values_list('key__name', 'value').iterator())

                result.append((info, col_label, name, criterion, l, gene_name, analysis_id))

            return result
        else:
            scaled_pvals = scale_df(list_enrichment_pvals)
            sns_heatmap = draw_heatmap(scaled_pvals, vmin=lower, vmax=upper)

            if save_file:
                sns_heatmap.savefig(pickledir + '/heatmap.svg')
                plt.close()

            return sns_heatmap

    return list_enrichment_pvals


def gene_list_enrichment_json(pickledir) -> Dict[str, Any]:
    df = gene_list_enrichment(
        pickledir,
        draw=False
    )
    names, criteria, analysis_ids, ls, uids = zip(*df.index)
    analyses = Analysis.objects.filter(pk__in=analysis_ids).prefetch_related('analysisdata_set',
                                                                             'analysisdata_set__key',
                                                                             'tf')

    result = []

    for (name, criterion, analysis_id, l, uid), *row in df.itertuples(name=None):
        info = {'filter': criterion, 'targets': l}
        try:
            analysis = analyses.get(pk=analysis_id)
            info['name'] = analysis.tf.gene_name_symbol
            info.update(analysis.meta_dict)
        except Analysis.DoesNotExist:
            pass

        result.append([info] + row)

    return {
        'columns': df.columns,
        'result': result
    }
