import os
from collections import defaultdict

import matplotlib

matplotlib.use('SVG')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import hypergeom


def heatmap(pickledir, cutoff, background):
    read_pickled_targetdbout(pickledir, cutoff, background)  # Func. to write


##################################################################
# function to read TargetDB output dataframe from pickles and
# convert into list of targets for each tF
def read_pickled_targetdbout(pickledir, cutoff=10, background=28775, save_file=True):
    pickled_pandas = pd.read_pickle(pickledir + '/' + 'tabular_output.pkl')
    anno_df = pd.read_pickle(pickledir + '/' + 'df_jsonanno.pkl')
    # type2_set is a nested dict, storing number of target genes for each exp and analysis
    # I had to store this data before filtering dataframe for target genes
    type2_set_dict= pd.read_pickle(pickledir+ '/' + 'df_eachtf_tgcount.pkl')

    # default cutoff is 10
    if not cutoff:
        cutoff=10
    if not background:
        background=28775
    # default background:
    # A. thaliana columbia tair10 genome(28775 genes(does not include transposable elements and pseudogenes))

    # raising exception here if target genes are not uploaded by the user
    try:
        pickled_targetgenes = pd.read_pickle(pickledir + '/' + 'df_targetgenes.pkl')
    except FileNotFoundError as e:
        raise FileNotFoundError('No target genes uploaded') from e

    # save targets for each TF into a dict
    targets_eachtf = defaultdict(list)

    # step2 ignore annotation and dap columns from the df
    col_ignore = ['Full Name', 'Family', 'Type', 'Name', 'List', 'UserList', 'TF Count']
    subset = pickled_pandas.iloc[:, ~(pickled_pandas.columns.get_level_values(0).isin(col_ignore))]

    expid_list = subset.columns.tolist()

    for val_expid in subset.columns.tolist():
        if 'Edges' in val_expid:
            target_eachanalysis = subset.loc[subset[(val_expid[0], val_expid[1], 'Edges')].notnull()]['ID'][
                                  3:].values
            genename = anno_df[val_expid[0].split('_')[0]][0]
            if genename == '-':
                genename = val_expid[0].split('_')[0]
            targets_eachtf[genename + ' || ' + val_expid[0] + ' || ' + val_expid[1]] = [l[0] for l in
                                                                                        target_eachanalysis]

    # Get data from modules
    module_names = list()

    # If a gene is in multiple lists, separate the lists by spaces and store all the names of user target gene lists
    # in a list
    # This comes with a restriction that user should not include any spaces in loaded targetgene lists.
    module_names_list = list(set(" ".join(pickled_targetgenes.List___UserList).split(" ")))
    # Replacing the dataframe spaces with $. If a gene is in multiple lists. These were separated by spaces
    pickled_targetgenes.replace(r'\s+', '$', regex=True, inplace=True)

    # empty pandas dataframe for overlaps and pvalues. Dataframes rows= analysis and columns= modules
    df_forheatmap = pd.DataFrame(np.nan, index=list(targets_eachtf.keys()), columns=module_names_list)
    dfpval_forheatmap = pd.DataFrame(np.nan, index=list(targets_eachtf.keys()), columns=module_names_list)

    rownamedict = dict()
    colnamedict = dict()

    # For loop for each user loaded list of target genes
    for val_module in module_names_list:
        # Get the genes in each module
        eachmodule_tg = pickled_targetgenes[(pickled_targetgenes['List___UserList'] == val_module) |
                                            (pickled_targetgenes['List___UserList'].str.startswith(val_module + '$')) |
                                            (pickled_targetgenes['List___UserList'].str.endswith('$' + val_module)) |
                                            (pickled_targetgenes['List___UserList'].str.contains(
                                                '\$' + val_module + '\$'))].index.tolist()
        # If the length of a loaded user list is less than 10 (by default), do not include this on the heatmap.
        # However this will appear on the tabular output
        if len(set(eachmodule_tg)) >= int(cutoff):
            colnamedict[val_module] = val_module + ' (' + str(len(set(eachmodule_tg))) + ')'

            for val_tg in targets_eachtf.keys():
                # If for an experiment, the number of target genes remaining are less than 10 (by default). By
                # remaining I mean
                #  when user upload a list of target genes then
                # Exclude this from
                # the heatmap. ## this is not very correct check the bitbucket issue #30
                if len(set(targets_eachtf[val_tg])) >= int(cutoff):

                    intersect_tg_mod = len(list(set(eachmodule_tg) & set(targets_eachtf[val_tg])))
                    df_forheatmap.ix[val_tg, val_module] = intersect_tg_mod  # assigning values to the dataframe
                    # The background here is from virtualplant
                    # Set one as a default and ask user for a parameter if they want to change the background
                    type_set2_length= type2_set_dict[val_tg.split('||')[1].strip()][val_tg.split('||')[2].strip()]
                    pval_uppertail = hypergeom.sf(intersect_tg_mod, int(background), len(set(eachmodule_tg)),
                                                  type_set2_length)
                    dfpval_forheatmap.ix[val_tg, val_module] = pval_uppertail
                    rownamedict[val_tg] = val_tg + ' (' + str(type_set2_length) + ')'

    # drops columns where all the values are nan. Modules with less than 10 genes
    dfpval_forheatmap.dropna(axis=1, how='all', inplace=True)
    # drops rows where all the values are nan. TFs with less than 10 targets for uploaded target list
    dfpval_forheatmap.dropna(axis=0, how='all', inplace=True)

    # wr2 = pd.ExcelWriter('raw_output_pval.xlsx')
    # df_forheatmap.to_excel(wr2, 'Sheet1')
    # wr2.save()

    dfpval_forheatmap.rename(columns=colnamedict, inplace=True)
    dfpval_forheatmap.rename(index=rownamedict, inplace=True)

    dfpval_forheatmap[dfpval_forheatmap < 1e-30] = 1e-30
    scaleddfpval_forhmap = -1 * np.log10(dfpval_forheatmap)
    scaleddfpval_forhmap.replace(np.inf, 1000, inplace=True)

    '''
    #hypergeom.sf(100, 12000, 3000, 400) is equal to 1-phyper(100,3000,12000-3000,400)
    '''

    sns_heatmap = sns.clustermap(scaleddfpval_forhmap, cmap="YlGnBu", cbar_kws={'label': 'Enrichment(-log10 p)'})
    plt.setp(sns_heatmap.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    plt.setp(sns_heatmap.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)

    if save_file:
        # get the absolute path of the directory
        outdirpath = os.path.abspath(pickledir)
        dirpath = '/'.join(outdirpath.split('/')[:-1])
        sns_heatmap.savefig(dirpath + '/' + pickledir.split('/')[-1].replace('_pickle', ''))
        plt.show()
        print('Generated= ', (dirpath + '/' + pickledir.split('/')[-1].replace('_pickle', '.svg')))
    else:
        return sns_heatmap
