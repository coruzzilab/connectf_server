import json
import os

import numpy as np
import pandas as pd

__all__ = ('create_cytoscape_data',)


def create_cytoscape_data(pickledir):
    tgdb_df, anno_df, db_tf, outdirpath = read_pickled_json(pickledir)  # Func. to write
    create_json(tgdb_df, anno_df, db_tf, outdirpath)

    return outdirpath


###################################################################
# Read pickled json file
def read_pickled_json(pickledir):
    tgdb_df = pd.read_pickle(pickledir + '/' + 'df_jsondata.pkl')
    anno_df = pd.read_pickle(pickledir + '/' + 'df_jsonanno.pkl')
    db_tf = pd.read_pickle(pickledir + '/' + 'df_jsondbtf.pkl')

    # create an output directory for downloadable zip file
    outdir = pickledir.replace('_pickle', '_json')
    if not os.path.exists(outdir):  # create output directory
        os.makedirs(outdir)

    # I am replacing here '',' ','    ' empty spaces with np.nan
    # cross check this with Zach: can this match with a string trapped inside space
    tgdb_df.replace(r'^\s*$', np.nan, regex=True, inplace=True)

    # wr= pd.ExcelWriter('json_pd_test.xlsx')
    # tgdb_df.to_excel(wr)
    # wr.close()

    # get the absolute path of the directory
    outdirpath = os.path.abspath(outdir)

    return tgdb_df, anno_df, db_tf, outdirpath


###################################################################
# function to create a json/object file for TFs
def create_json(query_res_df, mid_tfname, db_tf, output):
    # function 3 to create the raw data for creating json file for whole
    # genome TFs- Raw data means nodes, edges, size of the nodes, edge weight
    output_x = output + '/genome'
    edges_gm_rawdata, out_tf_genome, tf_genome_matrix_concat = get_json_root_edges_genome \
        (query_res_df, mid_tfname)
    geneid_genome = list(set(tf_genome_matrix_concat.columns.tolist() +
                    tf_genome_matrix_concat.index.tolist())) # For database TF matrix cols and rows (genome TF) are nodes

    create_json_object(geneid_genome, edges_gm_rawdata, out_tf_genome, mid_tfname, output_x)

    # function 4 to create the raw data for creating json file for
    # TFs-Targets- Raw data means nodes, edges, size of the nodes, edge weight
    output_x = output + '/targets'
    edges_tg_rawdata, out_tg_genome, tf_tg_matrix_concat = get_json_root_edges_dbase_TGs(query_res_df)
    geneid_TF_targets = list(set(tf_tg_matrix_concat.columns.tolist() + tf_tg_matrix_concat.index.tolist()))
    # query the database for geneid to genename mapping
    create_json_object(geneid_TF_targets, edges_tg_rawdata, out_tg_genome, mid_tfname, output_x)

    # function 2- get the data for creating json file for dbase
    # TFs data means nodes, edges, size of the nodes, edge weight
    output_view1 = output + '/dbase_view1'
    output_view2 = output + '/dbase_view2'
    output_view3 = output + '/dbase_view3'
    edges_db, out_tf_dbase, tf_subset_matrix_final = get_json_root_edges_dbase(query_res_df, db_tf)
    geneid_dbase = list(set(tf_subset_matrix_final.columns.tolist() + tf_subset_matrix_final.index.tolist()))
    db_mid_tfname = dict(
        (k, mid_tfname[k]) for k in geneid_dbase if k in mid_tfname)  # get a subset from full Ath annotation
    create_json_object(geneid_dbase, edges_db, out_tf_dbase, db_mid_tfname, output_view1)
    create_json_object(geneid_dbase, edges_db, out_tf_genome, db_mid_tfname, output_view2)
    create_json_object(geneid_dbase, edges_db, out_tg_genome, db_mid_tfname, output_view3)


###################################################################
# unstack the df to create edges links
def create_edge_links(rs_final_tfsubset_copy):
    # transpose the df to get the outdegree- in original df outdegree is columns
    transposed_df = rs_final_tfsubset_copy.transpose()
    # make the sif like format: exclude where tfs has no connection (0 in df)
    stacked_tmp_subset_df = pd.DataFrame(transposed_df[transposed_df > 0].stack().reset_index())
    stacked_tmp_subset_df.columns = ['source', 'target', 'val']  # assign columns to the df
    stacked_tmp_subset_df.dropna(axis=0, how='any')
    stacked_tmp_subset_df['id'] = stacked_tmp_subset_df['source'] + '_' + \
                                  stacked_tmp_subset_df['target'] + '_' + \
                                  stacked_tmp_subset_df['val']  # create ID for edges (TFID_TargetID)
    stacked_tmp_subset_df['color_induced'] = stacked_tmp_subset_df['val'].apply(
        lambda x: '#669900' if 'INDUCED' in x else '')
    stacked_tmp_subset_df['color_repressed'] = stacked_tmp_subset_df[
        'val'].apply(lambda x: '#FF6633' if 'REPRESSED' in x else '')
    stacked_tmp_subset_df['color_chipseq'] = stacked_tmp_subset_df['val'].apply(
        lambda x: '#666666' if x.isdigit() else '')
    stacked_tmp_subset_df['color'] = stacked_tmp_subset_df['color_chipseq'].map(
        str) + stacked_tmp_subset_df['color_induced'].map(str) + stacked_tmp_subset_df['color_repressed'].map(str)

    return stacked_tmp_subset_df


###################################################################
# Gets the root and edges for TFs present in the database
def get_json_root_edges_dbase(query_res_df, db_tf):
    # get rows with TF in the database and queried by user
    rs_final_tfsubset = query_res_df[query_res_df.index.isin(db_tf)]
    tf_subset_matrix = rs_final_tfsubset.notnull() * 1  # convert edges (values in df) into binary format
    tf_subset_matrix.rename(columns=lambda x: x.split('_')[0],
                            inplace=True)  # rename the df with TF name instead of expID
    tf_subset_matrix_final = tf_subset_matrix.groupby(tf_subset_matrix.columns, axis=1). \
        apply(lambda p: p.apply(lambda q: sum(q.values), axis=1))  # concat dfs for multiple TFs in df cols
    # Replacing here all the > 0 values with 1 to count no. of genes each TF is targeting
    tf_subset_matrix_final.where(tf_subset_matrix_final == 0, 1,
                                 inplace=True)  # Where keeps if is true and replaces false ones
    # Get edges data
    rs_final_tfsubset.rename(columns=lambda x: x.split('_')[0],
                             inplace=True)  # rename the df with TFname instead of expID
    # unstacking the df to create edges links
    stacked_tmp_subset_df = create_edge_links(rs_final_tfsubset)
    out_tf_dbase = (tf_subset_matrix_final * 1).sum(axis=0)  # store this data in a dict

    return stacked_tmp_subset_df, out_tf_dbase, tf_subset_matrix_final


######################################################################################
# Gets the node sizes for TFs based on number of TFs a TF is targeting in Ath genome
def get_json_root_edges_genome(rs_final_tfcol_subset, mid_tfname):
    list_tf_genome = [k for k in mid_tfname if
                      'TXNFACTOR' in mid_tfname[k]]  # get a subset from full Ath annotation
    # get the rows with TF (genome-wide) in database (means excluding targets that are not TFs)
    rs_tfgenome_sub = rs_final_tfcol_subset[rs_final_tfcol_subset.index.isin(list_tf_genome)]
    tf_genome_matrix = rs_tfgenome_sub.notnull() * 1  # convert edges (values in df) into binary format
    tf_genome_matrix.rename(columns=lambda x: x.split('_')[0], inplace=True)  # rename cols (TFs+expID) with TFs
    # merge cols for tfs (1 for interaction with row TF else 0)
    tf_genome_matrix_concat = tf_genome_matrix.groupby(tf_genome_matrix.columns, axis=1). \
        apply(lambda p: p.apply(lambda q: 1 if sum(q.values) > 0 else 0, axis=1))

    # Replacing here all the > 0 values with 1 to count no. of genes it is each TF is targeting
    tf_genome_matrix_concat.where(tf_genome_matrix_concat == 0, 1, inplace=True)
    # Get edges data
    rs_tfgenome_sub_copy = rs_tfgenome_sub.copy(deep=True)
    rs_tfgenome_sub_copy.rename(columns=lambda x: x.split('_')[0], inplace=True)  # rename the dataframe with TFname
    # unstacking the df to create edges links
    stacked_tmp_subset_df_gm = create_edge_links(rs_tfgenome_sub_copy)

    for stsd in tf_genome_matrix_concat:
        tf_genome_matrix_concat.loc[stsd, stsd] = 0  # replace values for TF-self connections to zero

    # Check point- It is converting all the values to binary. Be careful in
    # case sum as data is merged from multiple cols. So far I think it is fine.
    out_tf_genome = (tf_genome_matrix_concat * 1).sum(axis=0)  # store this data in a dict

    return stacked_tmp_subset_df_gm, out_tf_genome, tf_genome_matrix_concat


#####################################################################################
# Gets the root and edges for TFs and Targets present in the database
def get_json_root_edges_dbase_TGs(rs_final_tfcol_subset):
    tf_tg_matrix = rs_final_tfcol_subset.notnull() * 1  # convert edges (values in df) into binary format
    tf_tg_matrix.rename(columns=lambda x: x.split('_')[0], inplace=True)  # rename cols (TFs+experimentID) with TFs
    # merge cols for tfs (1 for interaction with row TF else 0)
    tf_tg_matrix_concat = tf_tg_matrix.groupby(tf_tg_matrix.columns, axis=1). \
        apply(lambda p: p.apply(lambda q: 1 if sum(q.values) > 0 else 0, axis=1))

    # Replacing here all the > 0 values with 1 to count no. of genes it is each TF is targeting
    tf_tg_matrix_concat.where(tf_tg_matrix_concat == 0, 1, inplace=True)
    # Get edges data
    rs_final_tfcol_copy = rs_final_tfcol_subset.copy(deep=True)
    # rename the dataframe with TF name instead of experimentIDs
    rs_final_tfcol_copy.rename(columns=lambda x: x.split('_')[0], inplace=True)
    # unstacking the df to create edges links
    stacked_tmp_subset_df_tg = create_edge_links(rs_final_tfcol_copy)

    edges_json_tg = list()  # object to be included in json final dict- that will be dupmed to a json file
    # converting df to a dict to create json format list
    for rowdf in (stacked_tmp_subset_df_tg.to_dict(orient='records')):
        edges_json_tg_dict = dict()
        edges_json_tg_dict['data'] = rowdf
        edges_json_tg.append(edges_json_tg_dict)

    out_tg_genome = (tf_tg_matrix_concat * 1).sum(axis=0)

    return stacked_tmp_subset_df_tg, out_tg_genome, tf_tg_matrix_concat


###################################################################
# create a json oject
# @profile
def create_json_object(geneid_x, edges_rawdata, out_tf_x, db_mid_tfname, output_x):
    elements = dict()  # CHANGE THIS TO A DICT WITH NODES AND EDGES
    edges_rawdata_copy = edges_rawdata.copy(deep=True)
    nodes = list()
    fontsize = out_tf_x.max()  # font size= size of TF with maximum no. of targets
    # font size edges= size of TF with maxno. of targets- normalize on log2 scale
    fontsize_edges = (np.log2(out_tf_x.max() + 2)) * 10

    # Assigning weights to the edges- this is only for chip-seq data
    edges_rawdata_copy['weight'] = edges_rawdata_copy['val'].apply(lambda x: (
        (x.count("1")) * fontsize_edges) if x.isdigit() else fontsize_edges)
    edges_rawdata_copy.drop(
        ['val', 'color_induced', 'color_repressed', 'color_chipseq'], 1,
        inplace=True)  # we don't need the dataframe values (0 1)

    edges_json_list = list()  # object to be included in json final dict-
    # that will be dupmed to a json file
    for rowdf in (
        edges_rawdata_copy.to_dict(orient='records')):  # converting df to a dict to create json format list
        edges_json_list_dict = dict()
        edges_json_list_dict['data'] = rowdf
        edges_json_list.append(edges_json_list_dict)

    for col_tfs in geneid_x:
        tmp_tfdict = dict()
        tmp_tfdict['data'] = dict()
        tmp_tfdict['data']['id'] = col_tfs  # gene ID

        if col_tfs in db_mid_tfname:
            if db_mid_tfname[col_tfs][1] == 'TXNFACTOR':
                tmp_tfdict['data']['type'] = db_mid_tfname[col_tfs][1]  # gene type
                tmp_tfdict['data']['color'] = "#00FF00"  # Assign color to TFs
                tmp_tfdict['data']['shape'] = "triangle"  # Assign shape to TFs
            elif db_mid_tfname[col_tfs][1] == 'PROTEIN_CODING':
                tmp_tfdict['data']['type'] = db_mid_tfname[col_tfs][1]  # gene type
                tmp_tfdict['data']['color'] = "#AED6F1"  # Assign color to non-TFs
                tmp_tfdict['data']['shape'] = "roundrectangle"  # Assign shape to non-TFs
            elif db_mid_tfname[col_tfs][1] == 'METABOLIC':
                tmp_tfdict['data']['type'] = db_mid_tfname[col_tfs][1]  # gene type
                tmp_tfdict['data']['color'] = "#D0ECE7"  # Assign color to non-TFs
                tmp_tfdict['data']['shape'] = "roundrectangle"  # Assign shape to non-TFs
            elif db_mid_tfname[col_tfs][1] == 'MOLECULE':
                tmp_tfdict['data']['type'] = db_mid_tfname[col_tfs][1]  # gene type
                tmp_tfdict['data']['color'] = "#FF9900"  # Assign color to non-TFs
                tmp_tfdict['data']['shape'] = "roundrectangle"  # Assign shape to non-TFs
            else:
                tmp_tfdict['data']['type'] = db_mid_tfname[col_tfs][1]  # gene type
                tmp_tfdict['data']['color'] = "#FF9900"  # Assign color to non-TFs
                tmp_tfdict['data']['shape'] = "roundrectangle"  # Assign shape to non-TFs

            '''
            elif (db_mid_tfname[col_tfs][1] == 'MOLECULE' or db_mid_tfname[col_tfs][1] == 'TRANSPOSABLE_ELEMENT_GENE' or
                  db_mid_tfname[col_tfs][1] == 'antisense_long_noncoding_rna' or
                  db_mid_tfname[col_tfs][1] == 'antisense_rna' or db_mid_tfname[col_tfs][1] == 'long_noncoding_rna' or
                  db_mid_tfname[col_tfs][1] == 'mirna' or db_mid_tfname[col_tfs][1] == 'novel_transcribed_region' or
                  db_mid_tfname[col_tfs][1] == 'other_rna' or db_mid_tfname[col_tfs][1] == 'pre_trna' or
                  db_mid_tfname[col_tfs][1] == 'pseudogene' or db_mid_tfname[col_tfs][1] == 'ribosomal_rna' or
                  db_mid_tfname[col_tfs][1] == 'ribosomal_rna' or db_mid_tfname[col_tfs][1] == 'small_nuclear_rna' or
                  db_mid_tfname[col_tfs][1] == 'small_nucleolar_rna'):
                tmp_tfdict['data']['type'] = db_mid_tfname[col_tfs][1]  # gene type
                tmp_tfdict['data']['color'] = "#FF9900"  # Assign color to non-TFs
                tmp_tfdict['data']['shape'] = "roundrectangle"  # Assign shape to non-TFs '''

        # assign node position
        tmp_tfdict['renderedPosition'] = dict()
        tmp_tfdict['renderedPosition']['x'] = 100
        tmp_tfdict['renderedPosition']['y'] = 100

        # assign node size (width and height)
        tmp_tfdict['style'] = dict()
        try:
            ht_wt = out_tf_x[col_tfs]
        except Exception:
            ht_wt = 0
        tmp_tfdict['style']['width'] = (np.log2(ht_wt + 2)) * 100
        tmp_tfdict['style']['height'] = (np.log2(ht_wt + 2)) * 100
        # for font size 1/3 of size of TF with max no. of targets
        tmp_tfdict['data']['weight'] = ((np.log2(fontsize + 2)) * 100) / 3

        if col_tfs in db_mid_tfname:
            if db_mid_tfname[col_tfs][0] == '-':
                tmp_tfdict['data']['name'] = col_tfs + ' (' + str(int(ht_wt)) + ')'  # gene name
            else:
                tmp_tfdict['data']['name'] = db_mid_tfname[col_tfs][0] + ' (' + str(int(ht_wt)) + ')'  # gene name
        else:
            tmp_tfdict['data']['name'] = col_tfs + ' (' + str(int(ht_wt)) + ')'  # gene name

        # assign group type
        tmp_tfdict['group'] = dict()
        tmp_tfdict['group'] = 'nodes'

        nodes.append(tmp_tfdict)  # append the style dict to a element list

    elements['nodes'] = nodes
    elements['edges'] = edges_json_list

    json_output_dict = dict()
    json_output_dict['elements'] = dict()  # create elements dict
    json_output_dict['elements'] = elements
    json_output_dict['layout'] = dict()  # Set layout
    json_output_dict['layout']['name'] = 'preset'

    if output_x:
        dir_path = os.path.dirname(os.path.realpath(output_x))
        with open(dir_path + '/' + output_x.split('/')[-1] + '_cy.json', 'w') as out_jsonfile:
            json.dump(json_output_dict, out_jsonfile, sort_keys=True, indent=4, ensure_ascii=False)
