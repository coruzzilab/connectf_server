import io
import os
from collections import defaultdict

import matplotlib

matplotlib.use('SVG')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import hypergeom


def heatmap(pickledir):
    read_pickled_targetdbout(pickledir)  # Func. to write


##################################################################
# function to read TargetDB output dataframe from pickles and
# convert into list of targets for each tF
def read_pickled_targetdbout(pickledir, save_file=True):
    pickled_pandas = pd.read_pickle(pickledir + '/' + 'tabular_output.pkl')
    anno_df = pd.read_pickle(pickledir + '/' + 'df_jsonanno.pkl')
    print('anno_df= ',type(anno_df))
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

    # expid_list= subset.columns.levels[0].tolist()

    expid_list = subset.columns.tolist()

    for val_expid in subset.columns.tolist():
        if 'Edges' in val_expid:
            target_eachanalysis = subset.loc[subset[(val_expid[0], val_expid[1], 'Edges')].notnull()]['ID'][
                                  3:].values
            genename= anno_df[val_expid[0].split('_')[0]][0]
            if genename == '-':
                genename = val_expid[0].split('_')[0]
            targets_eachtf[genename+'_'+val_expid[0]+'_'+val_expid[1]] = [l[0] for l in target_eachanalysis]

    # print('targets_eachtf= ',len(targets_eachtf.keys()))

    # Get data from modules
    module_names = list()
    # in case a target is in multiple lists, multiple lists will be separated by comma.
    # However for modules, this is not required as gene list is unique for different modules
    '''
    for x_tmp in pickled_targetgenes.List___UserList.unique():
        if ' ' in x_tmp.strip():
            for x_tmp_tmp in x_tmp.strip().split(' '):
                    module_names.append(x_tmp_tmp)
        else:
            module_names.append(x_tmp.strip())

    module_names_list= list(set(module_names))

    print('module_names_list= ',len(module_names_list))
    '''
    # Alternative: Testing is required for this part
    # If a gene is in multiple lists. Separate the space lists by
    module_names_list = list(set(" ".join(pickled_targetgenes.List___UserList).split(" ")))
    # print('module_names_list= ',len(module_names_list))
    pickled_targetgenes.replace(r'\s+', '$', regex=True, inplace=True)

    #wrs= pd.ExcelWriter('pickled_tggenes.xlsx')
    ##pickled_targetgenes.to_excel(wrs)
    #wrs.save()

    # print('pickled_targetgenes= ',pickled_targetgenes)

    # empty pandas dataframe with rows= analysis and columns= modules
    df_forheatmap = pd.DataFrame(np.nan, index=list(targets_eachtf.keys()), columns=module_names_list)
    dfpval_forheatmap = pd.DataFrame(np.nan, index=list(targets_eachtf.keys()), columns=module_names_list)
    # print('df_forheatmap.columns= ',df_forheatmap.columns)
    # print('df_forheatmap.index= ', df_forheatmap.index)

    rownamedict= dict()
    colnamedict= dict()

    for val_module in module_names_list:
        #print('val_module= ',val_module)
        eachmodule_tg = pickled_targetgenes[(pickled_targetgenes['List___UserList'] == val_module) |
                                            (pickled_targetgenes['List___UserList'].str.startswith(val_module + '$')) |
                                            (pickled_targetgenes['List___UserList'].str.endswith('$' + val_module)) |
                                            (pickled_targetgenes['List___UserList'].str.contains('\$' + val_module + '\$'))]. \
                                            index.tolist()

        colnamedict[val_module]= val_module+' ('+str(len(set(eachmodule_tg)))+')'

        for val_tg in targets_eachtf.keys():
            rownamedict[val_tg]= val_tg+' ('+str(len(set(targets_eachtf[val_tg])))+')'
            intersect_tg_mod = len(list(set(eachmodule_tg) & set(targets_eachtf[val_tg])))
            df_forheatmap.ix[val_tg, val_module] = intersect_tg_mod  # assigning values to the dataframe
            #print('intersect_tg_mod= ',intersect_tg_mod)
            #print('eachmodule_tg= ', len(set(eachmodule_tg)))
            #print('targets_eachtf[val_tg]= ',len(set(targets_eachtf[val_tg])))
            pval_uppertail = hypergeom.sf(intersect_tg_mod, 27655, len(set(eachmodule_tg)),
                                          len(set(targets_eachtf[val_tg])))
            dfpval_forheatmap.ix[val_tg, val_module] = pval_uppertail

    dfpval_forheatmap.rename(columns=colnamedict, inplace=True)
    dfpval_forheatmap.rename(index=rownamedict, inplace=True)

    #wr2 = pd.ExcelWriter('raw_output_pval.xlsx')
    #dfpval_forheatmap.to_excel(wr2, 'Sheet1')
    #wr2.save()

    dfpval_forheatmap[dfpval_forheatmap < 1e-30] = 1e-30
    scaleddfpval_forhmap = -1 * np.log10(dfpval_forheatmap)
    scaleddfpval_forhmap.replace(np.inf, 1000, inplace=True)

    if save_file:
        writer = pd.ExcelWriter('output_scaledpval.xlsx')
        scaleddfpval_forhmap.to_excel(writer, 'Sheet1')
        writer.save()

    '''
    writer = pd.ExcelWriter('output.xlsx')
    df_forheatmap.to_excel(writer, 'Sheet1')
    writer.save()

    writer = pd.ExcelWriter('output_pval.xlsx')
    dfpval_forheatmap.to_excel(writer, 'Sheet1')
    writer.save()

    #hypergeom.sf(100, 12000, 3000, 400) is equal to 1-phyper(100,3000,12000-3000,400)


    # create an empty dataframe for storing p-values
    # dataframe should be of dimension (number of analysis x number of modules)
    #df = pd.DataFrame(np.nan, index=[0, 1, 2, 3], columns=['A'])
    #df.ix['rowname', 'colname'] = 5.0
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
