from itertools import product

import matplotlib
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

matplotlib.use('SVG')
import matplotlib.pyplot as plt

import seaborn as sns

from queryapp.utils.parser import ANNOTATIONS


def scale_df(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df[df < 1e-30] = 1e-30
    df = -np.log10(df)
    df.replace(np.inf, 1000, inplace=True)
    return df


def draw_heatmap(df: pd.DataFrame):
    row_num, col_num = df.shape
    sns_heatmap = sns.clustermap(df,
                                 cmap="YlGnBu",
                                 cbar_kws={'label': 'Enrichment(-log10 p)'},
                                 col_cluster=col_num > 1,
                                 row_cluster=row_num > 1)
    plt.setp(sns_heatmap.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    plt.setp(sns_heatmap.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)

    return sns_heatmap


def save_heatmap(graph, pickledir):
    graph.savefig(pickledir + '/heatmap.svg')


##################################################################
# function to read TargetDB output dataframe from pickles and
# convert into list of targets for each tF
def heatmap(pickledir, cutoff=10, background=28775, draw=True, save_file=True):
    # raising exception here if target genes are not uploaded by the user
    try:
        name_to_list, list_to_name = pd.read_pickle(pickledir + '/target_genes.pickle.gz')
    except FileNotFoundError as e:
        raise FileNotFoundError('No target genes uploaded') from e

    query_result = pd.read_pickle(pickledir + '/tabular_output.pickle.gz')
    # type2_set is a nested dict, storing number of target genes for each exp and analysis
    # I had to store this data before filtering dataframe for target genes

    # default cutoff is 10
    # default background:
    # A. thaliana columbia tair10 genome(28775 genes(does not include transposable elements and pseudogenes))

    query_result = query_result.loc[:, (slice(None), slice(None), slice(None), 'EDGE')]

    # save targets for each TF into a dict
    targets = {}

    for val_expid, column in query_result.iteritems():
        target_eachanalysis = set(column.dropna().index)
        gene_id = val_expid[0]
        genename = ANNOTATIONS.at[gene_id, 'Name'] or gene_id
        targets['{0} || {1[1]} || {1[2]} ({2})'.format(
            genename, val_expid, len(target_eachanalysis))] = target_eachanalysis

    # @todo: wtf, clean this up
    # empty pandas dataframe for overlaps and pvalues. Dataframes rows= analysis and columns= modules
    dfpval_forheatmap = pd.DataFrame(np.nan, index=targets.keys(), columns=list_to_name.keys(), dtype=np.float64)

    colnames = {}

    for (analysis_name, analysis_list), (name, user_list) in product(targets.items(), list_to_name.items()):
        intersect_len = len(analysis_list & user_list)
        colnames[name] = "{} ({})".format(name, len(user_list))

        odds, pvalue = fisher_exact([
            [intersect_len, len(user_list - analysis_list)],
            [len(analysis_list - user_list), background - len(user_list | analysis_list)]
        ])
        dfpval_forheatmap.at[analysis_name, name] = pvalue

    dfpval_forheatmap.rename(columns=colnames, inplace=True)

    if draw:
        scaleddfpval_forhmap = scale_df(dfpval_forheatmap)
        sns_heatmap = draw_heatmap(scaleddfpval_forhmap)

        if save_file:
            save_heatmap(sns_heatmap, pickledir)
        else:
            return sns_heatmap

    else:
        return dfpval_forheatmap
