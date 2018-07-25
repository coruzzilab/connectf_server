from itertools import product
from typing import Optional

import matplotlib
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as hierarchy
from scipy.stats import fisher_exact

from queryapp.utils.parser import ANNOTATIONS

matplotlib.use('SVG')
import matplotlib.pyplot as plt

import seaborn as sns


def scale_df(df: pd.DataFrame, upper=30, lower=None) -> pd.DataFrame:
    df = df.copy()
    df = -np.log10(df)

    return df.clip(lower=lower, upper=upper)


def draw_heatmap(df: pd.DataFrame):
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
                                 **opts)
    plt.setp(sns_heatmap.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    plt.setp(sns_heatmap.ax_heatmap.xaxis.get_majorticklabels(), rotation=270)

    return sns_heatmap


##################################################################
# function to read TargetDB output dataframe from pickles and
# convert into list of targets for each tF
def heatmap(pickledir, background: Optional[int] = None, draw=True, save_file=False, upper=30, lower=None):
    # raising exception here if target genes are not uploaded by the user
    try:
        name_to_list, list_to_name = pd.read_pickle(pickledir + '/target_genes.pickle.gz')
    except FileNotFoundError as e:
        raise FileNotFoundError('No target genes uploaded') from e

    if background is None:
        background = ANNOTATIONS.shape[0]

    query_result = pd.read_pickle(pickledir + '/tabular_output_unfiltered.pickle.gz')

    # default cutoff is 10
    # default background:
    # A. thaliana columbia tair10 genome(28775 genes(does not include transposable elements and pseudogenes))

    query_result = query_result.loc[:, (slice(None), slice(None), slice(None), 'EDGE')]

    # save targets for each TF into a dict
    targets = {}

    for val_expid, column in query_result.iteritems():
        target_eachanalysis = set(column.dropna().index)
        gene_id = val_expid[0]
        try:
            # temporary fix for wonky names
            genename = ANNOTATIONS.at[gene_id, 'Name'] or gene_id
        except KeyError:
            genename = gene_id

        targets['{0} || {1[1]} || {1[2]} ({2})'.format(
            genename, val_expid, len(target_eachanalysis))] = target_eachanalysis

    dfpval_forheatmap = pd.DataFrame(index=targets.keys(), columns=list_to_name.keys(), dtype=np.float64)

    colnames = {}

    for (analysis_name, analysis_list), (name, user_list) in product(targets.items(), list_to_name.items()):
        intersect_len = len(analysis_list & user_list)
        colnames[name] = "{} ({})".format(name, len(user_list))

        odds, pvalue = fisher_exact([
            [intersect_len, len(user_list - analysis_list)],
            [len(analysis_list - user_list), background - len(user_list | analysis_list)]
        ], alternative='greater')
        dfpval_forheatmap.at[analysis_name, name] = pvalue

    dfpval_forheatmap.rename(columns=colnames, inplace=True)

    if draw:
        scaleddfpval_forhmap = scale_df(dfpval_forheatmap, upper=upper, lower=lower)
        sns_heatmap = draw_heatmap(scaleddfpval_forhmap)

        if save_file:
            sns_heatmap.savefig(pickledir + '/heatmap.svg')

        return sns_heatmap

    else:
        row_num, col_num = dfpval_forheatmap.shape
        if row_num > 1:
            z = hierarchy.linkage(dfpval_forheatmap.values, method='average', optimal_ordering=True)
            dfpval_forheatmap = dfpval_forheatmap.iloc[hierarchy.leaves_list(z), :]

        if col_num > 1:
            z = hierarchy.linkage(dfpval_forheatmap.values.T, method='average', optimal_ordering=True)
            dfpval_forheatmap = dfpval_forheatmap.iloc[:, hierarchy.leaves_list(z)]

        return dfpval_forheatmap
