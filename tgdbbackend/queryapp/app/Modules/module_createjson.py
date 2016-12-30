'''
This module returns three json objects:

	database TF- database TF
	database TF- whole genome TF
	database TF- whole genome (TF+nonTF)

Notes: create_json is the master function. Returns three objects.

Last updated: November 10, 2016
'''

##############
# Modules
import json
import os

import numpy as np
import pandas as pd
from create_mysqlDB import Genenames, Nodes


###################################################################
# function to create a json/object file for TFs
def create_json(sess, query_res_df, output):
    df_cols = query_res_df.columns.tolist()
    df_rows = query_res_df.index.tolist()

    rs_df_tf = sess.query(Nodes.node_name).filter(Nodes.node_type == 'TF'). \
        filter(Nodes.node_name.in_(
        df_rows)).all()  # get all the TFs in the database- target genes
    # should be in user query output dataframe
    df_tf_list = [x[0] for x in
                  rs_df_tf]  # convert database query output to list (TF list)

    # function 3 to create the raw data for creating json file for whole
    # genome TFs- Raw data means nodes, edges, size of the nodes, edge weight
    output_x = output + '_genome'
    edges_gm_rawdata, out_tf_genome, tf_genome_matrix_concat = \
        get_json_root_edges_genome(
        sess, query_res_df, df_rows)
    geneid_genome = list(set(
        tf_genome_matrix_concat.columns.tolist() +
        tf_genome_matrix_concat.index.tolist()))  # For database TF matrix I
    # need columns and rows (genome TF) as nodes
    mid_tfname_genome, tfname_tftype = id_name_mapping(sess,
                                                       geneid_genome)  #
    # query the database for geneid to genename mapping

    # json_object_WG= create_json_object(geneid_genome, edges_gm_rawdata,
    # out_tf_genome, mid_tfname_genome, tfname_tftype, output_x)

    # function 4 to create the raw data for creating json file for
    # TFs-Targets- Raw data means nodes, edges, size of the nodes, edge weight
    output_x = output + '_targets'
    edges_tg_rawdata, out_tg_genome, tf_tg_matrix_concat = \
        get_json_root_edges_dbase_TGs(
        query_res_df, df_rows)
    geneid_TF_targets = list(set(
        tf_tg_matrix_concat.columns.tolist() +
        tf_tg_matrix_concat.index.tolist()))
    mid_tfname_dbase, tfname_tftype = id_name_mapping(sess,
                                                      geneid_TF_targets)  #
    # query the database for geneid to genename mapping
    # json_object_TG= create_json_object(geneid_TF_targets, edges_tg_rawdata,
    #  out_tg_genome, mid_tfname_dbase, tfname_tftype, output_x)

    # function 2 to create the raw data for creating json file for database
    # TFs- Raw data means nodes, edges, size of the nodes, edge weight
    # output_view1= output+'_dbase_view1'
    output_view1 = output
    output_view2 = output + '_dbase_view2'
    output_view3 = output + '_dbase_view3'
    edge_rawdata, out_tf_dbase, tf_subset_matrix_final = \
        get_json_root_edges_dbase(
        query_res_df, df_tf_list)
    geneid_dbase = list(set(
        tf_subset_matrix_final.columns.tolist() +
        tf_subset_matrix_final.index.tolist()))  # For database TF matrix I
    # need both rows and columns as nodes
    mid_tfname_dbase, tfname_tftype = id_name_mapping(sess,
                                                      geneid_dbase)  # query
    # the database for geneid to genename mapping
    json_object_TFdbase_vw1 = create_json_object(geneid_dbase, edge_rawdata,
                                                 out_tf_dbase, mid_tfname_dbase,
                                                 tfname_tftype, output_view1)
    # json_object_TFdbase_vw2= create_json_object(geneid_dbase, edge_rawdata,
    #  out_tf_genome, mid_tfname_dbase, tfname_tftype, output_view2)
    # json_object_TFdbase_vw3= create_json_object(geneid_dbase, edge_rawdata,
    #  out_tg_genome, mid_tfname_dbase, tfname_tftype, output_view3)

    # function 1 to create the raw data for creating json file for query TFs-
    #  Raw data means nodes, edges, size of the nodes, edge weight
    output_query = output + '_query'
    edge_query_rawdata, out_query_dbase, tf_query_matrix_final = \
        get_json_root_edges_query(
        query_res_df, df_tf_list)
    geneid_query = list(set(
        tf_query_matrix_final.columns.tolist() +
        tf_query_matrix_final.index.tolist()))
    mid_tfname_query, tfname_tftype = id_name_mapping(sess,
                                                      geneid_query)  # query
    # the database for geneid to genename mapping
    # json_object_query= create_json_object(geneid_query, edge_query_rawdata,
    #  out_query_dbase, mid_tfname_query, tfname_tftype, output_query)

    # return json_object_WG, json_object_TG, json_object_query,
    # json_object_TFdbase_vw1, json_object_TFdbase_vw2, json_object_TFdbase_vw3
    return json_object_TFdbase_vw1


###################################################################
# unstack the df to create edges links
def create_edge_links(rs_final_tfsubset_copy):
    transposed_df = rs_final_tfsubset_copy.transpose()  # transpose the df to
    #  get the outdegree- because in original df outdegree is columns
    stacked_tmp_subset_df = pd.DataFrame(transposed_df[
                                             transposed_df > 0].stack(

    ).reset_index())  # make the sif like format: exclude where tfs has no
    # connection (0 in df)
    stacked_tmp_subset_df.columns = ['source', 'target',
                                     'val']  # assign columns to the df
    stacked_tmp_subset_df.dropna(axis=0, how='any')
    stacked_tmp_subset_df['id'] = stacked_tmp_subset_df['source'] + '_' + \
                                  stacked_tmp_subset_df['target'] + '_' + \
                                  stacked_tmp_subset_df[
                                      'val']  # create ID for edges (
    # TFID_TargetID)
    stacked_tmp_subset_df['color_induced'] = stacked_tmp_subset_df['val'].apply(
        lambda x: '#669900' if 'INDUCED' in x else '')
    stacked_tmp_subset_df['color_repressed'] = stacked_tmp_subset_df[
        'val'].apply(lambda x: '#FF6633' if 'REPRESSED' in x else '')
    stacked_tmp_subset_df['color_chipseq'] = stacked_tmp_subset_df['val'].apply(
        lambda x: '#666666' if x.isdigit() else '')
    stacked_tmp_subset_df['color'] = stacked_tmp_subset_df['color_chipseq'].map(
        str) + stacked_tmp_subset_df['color_induced'].map(str) + \
                                     stacked_tmp_subset_df[
                                         'color_repressed'].map(str)

    return stacked_tmp_subset_df


###################################################################
# Gets the root and edges for TFs present in the database
def get_json_root_edges_query(rs_final_tfcol_subset, df_tf_list):
    df_columns = rs_final_tfcol_subset.columns.tolist()
    cols_tfs_list = list(
        set([x.split('_')[0] for x in df_columns]))  # from df columns create a

    rs_final_tfquery = rs_final_tfcol_subset[rs_final_tfcol_subset.index.isin(
        cols_tfs_list)]  # get rows with TF queried by user
    tf_query_matrix = rs_final_tfquery.notnull() * 1  # convert edges (values
    #  in df) into binary format
    tf_query_matrix.rename(columns=lambda x: x.split('_')[0],
                           inplace=True)  # rename the dataframe with TF name
    #  instead of experimentID
    tf_query_matrix_final = tf_query_matrix.groupby(tf_query_matrix.columns,
                                                    axis=1). \
        apply(lambda p: p.apply(lambda q: sum(q.values),
                                axis=1))  # concat dataframes for multiple
    # TFs in dataframe columns
    # Replacing here all the > 0 values with 1 to count no. of genes each TF
    # is targeting
    tf_query_matrix_final.where(tf_query_matrix_final == 0, 1,
                                inplace=True)  # Where command keeps values
    # where condition is true and replaces false ones

    # Get edges data
    rs_final_tfquery_copy = rs_final_tfquery.copy(deep=True)
    rs_final_tfquery_copy.rename(columns=lambda x: x.split('_')[0],
                                 inplace=True)  # rename the dataframe with
    # TF name instead of experimentID
    # unstacking the df to create edges links
    stacked_query_subset_df = create_edge_links(rs_final_tfquery_copy)

    # identify the root node: discard TF's connection to itself
    # Check point- Here the tf_subset_matrix_final will be changed- TF to
    # itself connections will be removed to identify root node
    tf_query_matrix_final_copy = tf_query_matrix_final.copy(
        deep=True)  # creating a deep copy of the df to avoid changing the df

    for stsd in tf_query_matrix_final_copy:
        tf_query_matrix_final_copy.loc[
            stsd, stsd] = 0  # replace values for TF-self connections to zero

    out_querytf_dbase = (tf_query_matrix_final_copy * 1).sum(
        axis=0)  # store this data in a dict
    tf_query_matrix_final_copy['row_sum'] = tf_query_matrix_final_copy.sum(
        axis=1)
    tf_query_matrix_final_copy['column_sum'] = tf_query_matrix_final_copy.sum(
        axis=0)

    return stacked_query_subset_df, out_querytf_dbase, tf_query_matrix_final


###################################################################
# Gets the root and edges for TFs present in the database
def get_json_root_edges_dbase(rs_final_tfcol_subset, df_tf_list):
    rs_final_tfsubset = rs_final_tfcol_subset[rs_final_tfcol_subset.index.isin(
        df_tf_list)]  # get rows with TF in the database and queried by user
    tf_subset_matrix = rs_final_tfsubset.notnull() * 1  # convert edges (
    # values in df) into binary format
    tf_subset_matrix.rename(columns=lambda x: x.split('_')[0],
                            inplace=True)  # rename the dataframe with TF
    # name instead of experimentID
    tf_subset_matrix_final = tf_subset_matrix.groupby(tf_subset_matrix.columns,
                                                      axis=1). \
        apply(lambda p: p.apply(lambda q: sum(q.values),
                                axis=1))  # concat dataframes for multiple
    # TFs in dataframe columns
    uniqueele = list(set(tf_subset_matrix.columns.tolist()) - set(
        tf_subset_matrix.index.tolist()))  # get gene id of tfs that are
    # present in df columns but not rows
    # Replacing here all the > 0 values with 1 to count no. of genes each TF
    # is targeting
    tf_subset_matrix_final.where(tf_subset_matrix_final == 0, 1,
                                 inplace=True)  # Where command keeps values
    # where condition is true and replaces false ones

    # Get edges data
    rs_final_tfsubset_copy = rs_final_tfsubset.copy(deep=True)
    rs_final_tfsubset_copy.rename(columns=lambda x: x.split('_')[0],
                                  inplace=True)  # rename the dataframe with
    # TF name instead of experimentID
    # unstacking the df to create edges links
    stacked_tmp_subset_df = create_edge_links(rs_final_tfsubset_copy)

    # identify the root node: discard TF's connection to itself
    # Check point- Here the tf_subset_matrix_final will be changed- TF to
    # itself connections will be removed to identify root node
    tf_subset_matrix_final_copy = tf_subset_matrix_final.copy(
        deep=True)  # creating a deep copy of the df to avoid changing the df

    for stsd in tf_subset_matrix_final_copy:
        tf_subset_matrix_final_copy.loc[
            stsd, stsd] = 0  # replace values for TF-self connections to zero

    out_tf_dbase = (tf_subset_matrix_final_copy * 1).sum(
        axis=0)  # store this data in a dict
    tf_subset_matrix_final_copy['row_sum'] = tf_subset_matrix_final_copy.sum(
        axis=1)
    tf_subset_matrix_final_copy['column_sum'] = tf_subset_matrix_final_copy.sum(
        axis=0)

    return stacked_tmp_subset_df, out_tf_dbase, tf_subset_matrix_final


#####################################################################################
# Gets the node sizes for TFs based on number of TFs a TF is targeting in Ath
#  genome
def get_json_root_edges_genome(sess, rs_final_tfcol_subset, df_rows):
    rs_df_tf_genome = sess.query(Genenames.ath_id).filter(
        Genenames.ath_gene_type == 'TXNFACTOR'). \
        filter(Genenames.ath_id.in_(
        df_rows)).all()  # get genes from df rows that are defined as TFs Ath
    #  annotation- ath_gene_type in the database
    df_tf_genome = [x[0] for x in
                    rs_df_tf_genome]  # convert database query (genome)
    # output to list (TF list)
    rs_tfgenome_sub = rs_final_tfcol_subset[rs_final_tfcol_subset.index.isin(
        df_tf_genome)]  # get the rows with TF (genome-wide) in database (
    # means excluding targets that are not TFs)
    tf_genome_matrix = rs_tfgenome_sub.notnull() * 1  # convert edges (values
    #  in df) into binary format
    tf_genome_matrix.rename(columns=lambda x: x.split('_')[0],
                            inplace=True)  # rename cols (TFs+experimentID)
    # with TFs
    tf_genome_matrix_concat = tf_genome_matrix.groupby(tf_genome_matrix.columns,
                                                       axis=1). \
        apply(lambda p: p.apply(lambda q: 1 if sum(q.values) > 0 else 0,
                                axis=1))  # merge cols for tfs (1 for
    # interaction with row TF else 0)

    # Replacing here all the > 0 values with 1 to count no. of genes it is
    # each TF is targeting
    tf_genome_matrix_concat.where(tf_genome_matrix_concat == 0, 1,
                                  inplace=True)  # Where keeps values where
    # condition is true and replaces false ones
    # Get edges data
    rs_tfgenome_sub_copy = rs_tfgenome_sub.copy(deep=True)
    rs_tfgenome_sub_copy.rename(columns=lambda x: x.split('_')[0],
                                inplace=True)  # rename the dataframe with TF
    #  name instead of experimentID
    # unstacking the df to create edges links
    stacked_tmp_subset_df_gm = create_edge_links(rs_tfgenome_sub_copy)

    for stsd in tf_genome_matrix_concat:
        tf_genome_matrix_concat.loc[
            stsd, stsd] = 0  # replace values for TF-self connections to zero

    # Check point- It is converting all the values to binary. Be careful in
    # case sum as data is merged from multiple cols. So far I think it is fine.
    out_tf_genome = (tf_genome_matrix_concat * 1).sum(
        axis=0)  # store this data in a dict

    return stacked_tmp_subset_df_gm, out_tf_genome, tf_genome_matrix_concat


#####################################################################################
# Gets the root and edges for TFs and Targets present in the database
def get_json_root_edges_dbase_TGs(rs_final_tfcol_subset, df_tf_list):
    tf_tg_matrix = rs_final_tfcol_subset.notnull() * 1  # convert edges (
    # values in df) into binary format
    tf_tg_matrix.rename(columns=lambda x: x.split('_')[0],
                        inplace=True)  # rename cols (TFs+experimentID) with TFs
    tf_tg_matrix_concat = tf_tg_matrix.groupby(tf_tg_matrix.columns, axis=1). \
        apply(lambda p: p.apply(lambda q: 1 if sum(q.values) > 0 else 0,
                                axis=1))  # merge cols for tfs (1 for
    # interaction with row TF else 0)

    # Replacing here all the > 0 values with 1 to count no. of genes it is
    # each TF is targeting
    tf_tg_matrix_concat.where(tf_tg_matrix_concat == 0, 1,
                              inplace=True)  # Where keeps values where
    # condition is true and replaces false ones
    # Get edges data
    rs_final_tfcol_copy = rs_final_tfcol_subset.copy(deep=True)
    rs_final_tfcol_copy.rename(columns=lambda x: x.split('_')[0],
                               inplace=True)  # rename the dataframe with TF
    # name instead of experimentIDs
    # unstacking the df to create edges links
    stacked_tmp_subset_df_tg = create_edge_links(rs_final_tfcol_copy)

    edges_json_tg = list()  # object to be included in json final dict- that
    # will be dupmed to a json file
    for rowdf in (stacked_tmp_subset_df_tg.to_dict(
        orient='records')):  # converting df to a dict to create json format
        # list
        edges_json_tg_dict = dict()
        edges_json_tg_dict['data'] = rowdf
        edges_json_tg.append(edges_json_tg_dict)

    out_tg_genome = (tf_tg_matrix_concat * 1).sum(
        axis=0)  # store this data in a dict

    return stacked_tmp_subset_df_tg, out_tg_genome, tf_tg_matrix_concat


###################################################################
# create a json oject
def create_json_object(geneid_x, edges_rawdata, out_tf_x, mid_tfname,
                       tfname_tftype, output_x):
    elements = dict()  # CHANGE THIS TO A DICT WITH NODES AND EDGES
    # create a css dict for style sheet
    css_node = dict()
    css_node['content'] = 'data(name)'
    css_node['font-family'] = 'helvetica'
    css_node['font-size'] = 100
    css_node['font-weight'] = 'bold'
    css_node['text-valign'] = 'center'
    css_node['color'] = '#000000'
    css_node['shape'] = 'triangle'
    css_node['background-color'] = '#65CA7C'
    css_edge = dict()
    css_edge['width'] = 3
    css_edge['color'] = '#000000'
    css_edge['target-arrow-shape'] = 'triangle'
    css_edge['background-color'] = '#000000'
    edges_rawdata_copy = edges_rawdata.copy(deep=True)
    nodes = list()
    fontsize = out_tf_x.max()  # font size= size of TF with maximum no. of
    # targets
    fontsize_edges = (np.log2(
        out_tf_x.max() + 2)) * 10  # font size edges= size of TF with maximum
    #  no. of targets- normalize on log2 scale

    # Assigning weigths to the edges- this is only for chip-seq data
    edges_rawdata_copy['weight'] = edges_rawdata_copy['val'].apply(lambda x: (
    (x.count("1")) * fontsize_edges) if x.isdigit() else fontsize_edges)
    edges_rawdata_copy.drop(
        ['val', 'color_induced', 'color_repressed', 'color_chipseq'], 1,
        inplace=True)  # we don't need the dataframe values (0 1)

    edges_json_list = list()  # object to be included in json final dict-
    # that will be dupmed to a json file
    for rowdf in (edges_rawdata_copy.to_dict(
        orient='records')):  # converting df to a dict to create json format
        # list
        edges_json_list_dict = dict()
        edges_json_list_dict['data'] = rowdf
        edges_json_list.append(edges_json_list_dict)

    for col_tfs in geneid_x:
        tmp_tfdict = dict()
        tmp_tfdict['data'] = dict()
        tmp_tfdict['data']['id'] = col_tfs  # gene ID

        if col_tfs in tfname_tftype:
            if tfname_tftype[col_tfs] == 'TXNFACTOR':
                tmp_tfdict['data']['type'] = tfname_tftype[col_tfs]  # gene type
                tmp_tfdict['data']['color'] = "#00FF00"  # Assign color to TFs
                tmp_tfdict['data']['shape'] = "triangle"  # Assign shape to TFs
            elif tfname_tftype[col_tfs] == 'PROTEIN_CODING':
                tmp_tfdict['data']['type'] = tfname_tftype[col_tfs]  # gene type
                tmp_tfdict['data'][
                    'color'] = "#AED6F1"  # Assign color to non-TFs
                tmp_tfdict['data'][
                    'shape'] = "roundrectangle"  # Assign shape to non-TFs
            elif tfname_tftype[col_tfs] == 'METABOLIC':
                tmp_tfdict['data']['type'] = tfname_tftype[col_tfs]  # gene type
                tmp_tfdict['data'][
                    'color'] = "#D0ECE7"  # Assign color to non-TFs
                tmp_tfdict['data'][
                    'shape'] = "roundrectangle"  # Assign shape to non-TFs
            elif tfname_tftype[col_tfs] == 'MOLECULE':
                tmp_tfdict['data']['type'] = tfname_tftype[col_tfs]  # gene type
                tmp_tfdict['data'][
                    'color'] = "#FF9900"  # Assign color to non-TFs
                tmp_tfdict['data'][
                    'shape'] = "roundrectangle"  # Assign shape to non-TFs

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
        tmp_tfdict['data']['weight'] = ((np.log2(
            fontsize + 2)) * 100) / 3  # for font size 1/3 of size of TF with maximum no. of targets

        if col_tfs in mid_tfname:
            if mid_tfname[col_tfs] == '-':
                tmp_tfdict['data']['name'] = col_tfs + ' (' + str(
                    int(ht_wt)) + ')'  # gene name
            else:
                tmp_tfdict['data']['name'] = mid_tfname[col_tfs] + ' (' + str(
                    int(ht_wt)) + ')'  # gene name
        else:
            tmp_tfdict['data']['name'] = col_tfs + ' (' + str(
                int(ht_wt)) + ')'  # gene name

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
    json_object = json.dumps(json_output_dict)

    if output_x:
        dir_path = os.path.dirname(os.path.realpath(output_x))
        # mydir= '/Users/Reetu/Documents/Projects/TargetDB/tgdbbackend/tgdbbackend/static/queryBuilder'
        with open(dir_path + '/' + output_x.split('/')[-1] + '_cy.json',
                  'wb') as \
            out_jsonfile:
            # print '***= ','/Users/Reetu/Documents/Projects/TargetDB/tgdbbackend/tgdbbackend/static/queryBuilder'+'/'+output_x.split('/')[-1]+'_cy.json'
            json.dump(json_output_dict, out_jsonfile, sort_keys=True, indent=4,
                      ensure_ascii=False)

    return json_object


###################################################################
# Function to create gene_id to gene_name mapping
def id_name_mapping(sess, geneidlist):
    mid_tfname = dict()
    tfname_tftype = dict()
    for i_mid in geneidlist:

        tf_name = sess.query(Genenames.ath_id, Genenames.ath_name,
                             Genenames.ath_gene_type). \
            filter(Genenames.ath_id == i_mid).all()
        if tf_name:
            mid_tfname[tf_name[0][0]] = str(tf_name[0][1])
            tfname_tftype[tf_name[0][0]] = str(tf_name[0][2])
        else:
            mid_tfname[i_mid] = '-'
            tfname_tftype[i_mid] = 'unknown'

    return mid_tfname, tfname_tftype
