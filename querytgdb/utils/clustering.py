from collections import defaultdict

import matplotlib

matplotlib.use('SVG')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import fisher_exact


def heatmap(pickledir, cutoff, background):
    read_pickled_targetdbout(pickledir, cutoff, background)  # Func. to write


def scale_df(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df[df < 1e-30] = 1e-30
    df = -np.log10(df)
    df.replace(np.inf, 1000, inplace=True)
    return df


def draw_heatmap(df: pd.DataFrame):
    sns_heatmap = sns.clustermap(df, cmap="YlGnBu", cbar_kws={'label': 'Enrichment(-log10 p)'})
    plt.setp(sns_heatmap.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    plt.setp(sns_heatmap.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)

    return sns_heatmap


def save_heatmap(graph, pickledir):
    graph.savefig(pickledir + '/heatmap.svg')


##################################################################
# function to read TargetDB output dataframe from pickles and
# convert into list of targets for each tF
def read_pickled_targetdbout(pickledir, cutoff=10, background=28775, draw=True, save_file=True):
    # raising exception here if target genes are not uploaded by the user
    try:
        pickled_targetgenes = pd.read_pickle(pickledir + '/df_targetgenes.pkl')
    except FileNotFoundError as e:
        raise FileNotFoundError('No target genes uploaded') from e

    pickled_pandas = pd.read_pickle(pickledir + '/tabular_output.pkl')
    anno_df = pd.read_pickle(pickledir + '/df_jsonanno.pkl')
    # type2_set is a nested dict, storing number of target genes for each exp and analysis
    # I had to store this data before filtering dataframe for target genes
    type2_set_dict = pd.read_pickle(pickledir + '/df_eachtf_tgcount.pkl')

    # default cutoff is 10
    # default background:
    # A. thaliana columbia tair10 genome(28775 genes(does not include transposable elements and pseudogenes))

    # save targets for each TF into a dict
    targets_eachtf = defaultdict(list)

    # step2 ignore annotation and dap columns from the df
    col_ignore = ['Full Name', 'Family', 'Type', 'Name', 'List', 'UserList', 'TF Count']
    # print(pickled_pandas.loc[:, pickled_pandas.columns.get_level_values(0).isin(col_ignore)])
    subset = pickled_pandas.loc[:, ~pickled_pandas.columns.get_level_values(0).isin(col_ignore)]

    for val_expid in subset.columns[subset.columns.get_level_values(2) == 'Edges']:
        target_eachanalysis = subset.loc[subset[val_expid].notnull(), 'ID'].iloc[3:, 0].tolist()
        gene_id = val_expid[0].partition('_')[0]
        genename = anno_df.at[gene_id, 'Full Name___Gene Full Name']
        if genename == '-':
            genename = gene_id
        targets_eachtf['{0} || {1[0]} || {1[1]}'.format(genename, val_expid)] = target_eachanalysis

    # If a gene is in multiple lists, separate the lists by spaces and store all the names of user target gene lists
    # in a list
    # This comes with a restriction that user should not include any spaces in loaded targetgene lists.
    # @todo: wtf, clean this up
    module_names_list = np.unique(pickled_targetgenes.List___UserList)
    # Replacing the dataframe spaces with $. If a gene is in multiple lists. These were separated by spaces

    # empty pandas dataframe for overlaps and pvalues. Dataframes rows= analysis and columns= modules
    df_forheatmap = pd.DataFrame(np.nan, index=list(targets_eachtf.keys()), columns=module_names_list)
    dfpval_forheatmap = df_forheatmap.copy(deep=True)

    rownamedict = dict()
    colnamedict = dict()

    # For loop for each user loaded list of target genes
    for val_module in module_names_list:
        # Get the genes in each module
        eachmodule_tg = pickled_targetgenes[
            pickled_targetgenes['List___UserList'].str.contains(r'\b' + val_module + r'\b')].index.tolist()
        # If the length of a loaded user list is less than 10 (by default), do not include this on the heatmap.
        # However this will appear on the tabular output
        eachmodule_tg_len = len(set(eachmodule_tg))
        if eachmodule_tg_len >= cutoff:
            colnamedict[val_module] = f'{val_module} ({eachmodule_tg_len})'

            for val_tg, genes in targets_eachtf.items():
                # If for an experiment, the number of target genes remaining are less than 10 (by default). By
                # remaining I mean
                #  when user upload a list of target genes then
                # Exclude this from the heatmap.
                gene_name, experiment, analysis = val_tg.split(' || ')
                type_set2_length = type2_set_dict[experiment][analysis]
                if type_set2_length >= cutoff:
                    intersect_tg_mod = len(set(eachmodule_tg) & set(genes))
                    df_forheatmap.at[val_tg, val_module] = intersect_tg_mod  # assigning values to the dataframe
                    # The background here is from virtualplant
                    # Set one as a default and ask user for a parameter if they want to change the background
                    odds, pval_uppertail = fisher_exact([[intersect_tg_mod, eachmodule_tg_len - intersect_tg_mod],
                                                         [type_set2_length - intersect_tg_mod,
                                                          background - eachmodule_tg_len - type_set2_length +
                                                          intersect_tg_mod]],
                                                        alternative='greater')
                    dfpval_forheatmap.at[val_tg, val_module] = pval_uppertail
                    rownamedict[val_tg] = f'{val_tg} ({type_set2_length})'

    # drops columns where all the values are nan. Modules with less than 10 genes
    dfpval_forheatmap.dropna(axis=1, how='all', inplace=True)
    # drops rows where all the values are nan. TFs with less than 10 targets for uploaded target list
    dfpval_forheatmap.dropna(axis=0, how='all', inplace=True)

    dfpval_forheatmap.rename(columns=colnamedict, inplace=True)
    dfpval_forheatmap.rename(index=rownamedict, inplace=True)

    if draw:
        scaleddfpval_forhmap = scale_df(dfpval_forheatmap)
        sns_heatmap = draw_heatmap(scaleddfpval_forhmap)

        if save_file:
            save_heatmap(sns_heatmap, pickledir)
        else:
            return sns_heatmap

    else:
        return dfpval_forheatmap
