import sys
from itertools import count, product
from typing import Optional
from uuid import uuid4

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as hierarchy
import seaborn as sns
from scipy.stats import fisher_exact

from ..models import Analysis
from ..utils import NAME_REGEX, column_string
from ..utils.parser import ANNOTATIONS


def scale_df(df: pd.DataFrame) -> pd.DataFrame:
    df = -np.log10(df)
    return df.clip_upper(sys.maxsize)


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
                                 row_cluster=row_num > 1,
                                 col_cluster=col_num > 1,
                                 **kwargs,
                                 **opts)
    plt.setp(sns_heatmap.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    # plt.setp(sns_heatmap.ax_heatmap.yaxis.get_label(), text="Analyses", rotation=270)
    plt.setp(sns_heatmap.ax_heatmap.xaxis.get_majorticklabels(), rotation=270)
    plt.setp(sns_heatmap.ax_heatmap.xaxis.get_label(), text="Target Genes")

    return sns_heatmap


def heatmap(pickledir, background: Optional[int] = None, draw=True, legend=False, save_file=False, upper=None,
            lower=None):
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

    query_result = query_result.loc[:, (slice(None), slice(None), slice(None), 'EDGE')]
    query_result.columns = query_result.columns.droplevel(3)

    names, exp_ids, analysis_ids = zip(*query_result.columns)
    analyses = Analysis.objects.filter(
        name__in=analysis_ids,
        experiment__name__in=exp_ids
    ).prefetch_related('experiment')

    # save targets for each TF into a dict
    targets = {}

    for (name, exp_id, analysis_id), column in query_result.iteritems():
        analysis_targets = set(column.dropna().index)

        m = NAME_REGEX.match(name)

        if m:
            name, criterion = m.groups('')
        else:
            criterion = ''

        targets[(name, criterion, exp_id, analysis_id, len(analysis_targets), uuid4())] = analysis_targets

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

        for (name, criterion, exp_id, analysis_id, l, uid), col_name in orig_index:
            tf = analyses.get(name=analysis_id, experiment__name=exp_id).experiment.tf
            indices.append('{2} — {0}{1} ({3})'.format(tf.gene_id, ' ' + tf.name if tf.name else '', col_name, l))

        list_enrichment_pvals.index = indices

        if legend:
            result = []

            for idx, col_label in orig_index:
                name, criterion, exp_id, analysis_id, l, uid = idx

                info = {'name': name, 'filter': criterion}
                analysis = analyses.get(
                    name=analysis_id,
                    experiment__name=exp_id)

                gene_id = analysis.experiment.tf.gene_id

                try:
                    gene_name = ANNOTATIONS.at[gene_id, 'Name']
                except KeyError:
                    gene_name = ''

                info.update(analysis.analysisdata_set.values_list('key', 'value').iterator())
                info.update(analysis.experiment.experimentdata_set.values_list('key', 'value').iterator())

                result.append((info, col_label, name, criterion, l, gene_name, analysis_id))

            return result
        else:
            scaled_pvals = scale_df(list_enrichment_pvals)
            sns_heatmap = draw_heatmap(scaled_pvals, vmin=lower, vmax=upper)

            if save_file:
                sns_heatmap.savefig(pickledir + '/heatmap.svg')

            return sns_heatmap

    else:
        row_num, col_num = list_enrichment_pvals.shape
        if row_num > 1:
            z = hierarchy.linkage(list_enrichment_pvals.values, method='average', optimal_ordering=True)
            list_enrichment_pvals = list_enrichment_pvals.iloc[hierarchy.leaves_list(z), :]

        if col_num > 1:
            z = hierarchy.linkage(list_enrichment_pvals.values.T, method='average', optimal_ordering=True)
            list_enrichment_pvals = list_enrichment_pvals.iloc[:, hierarchy.leaves_list(z)]

        return list_enrichment_pvals
