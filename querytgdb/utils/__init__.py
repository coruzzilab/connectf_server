import os
import pickle
import random
import re
import string
import warnings
from collections import OrderedDict, defaultdict
from itertools import groupby, tee
from operator import itemgetter, methodcaller
from string import whitespace
from typing import Dict, Generator, Tuple

import numpy as np
import pandas as pd
from pandas.core.common import SettingWithCopyWarning
from pandas.errors import PerformanceWarning

from ..models import Analysis, AnalysisIddata, Annotation, Interactions, MetaIddata, Metadata, ReferenceId, \
    Regulation, TargetDBTF

__all__ = ('query_tgdb', 'rand_string')

warnings.simplefilter(action='ignore', category=(FutureWarning, PerformanceWarning, SettingWithCopyWarning))


def rand_string(length):
    return ''.join(random.choices(string.ascii_letters + string.digits + '_', k=length))


def get_p_values(orig_meta: str) -> Generator[float, None, None]:
    p_values = re.findall(r'p-value\s*=\s*(.+?)(?=\s|])', orig_meta, flags=re.I)

    for p in p_values:
        try:
            yield float(p)
        except ValueError:
            pass


def query_tgdb(tf_query, edges, metadata, target_genes, output):
    # check if the command line arguments provided are ok
    if tf_query is None or not any(map((lambda x: str.strip(x, '[]' + whitespace)), tf_query)):
        tf_query = ['[OR', '[ALLTF]]']
        # raise TypeError('Either generate a table for all the TFs (--t= OR [ALLTF] \n or '
        #                 'query TFs, both can not be none \n')

    p_value = None
    if metadata:
        orig_meta = ' '.join(metadata)

        try:
            p_value = min(get_p_values(orig_meta))
        except ValueError:
            pass

        # remove p-value for now
        metadata = re.sub(r'\b\s*(?:And|Or|AndNot)?\s*p-value\s*=.+?(?=\s|])', '', orig_meta, flags=re.I)
        metadata = re.sub(r'(?<=\[)\s*(?:And|Or|AndNot)\s*', '', metadata, flags=re.I)
        metadata = re.sub(r'\s*\[\s*\]', '', metadata, flags=re.I)
        while metadata:
            _metadata = re.sub(r'\s*\[\s*\]', '', metadata, flags=re.I)
            if _metadata == metadata:
                break
            else:
                metadata = _metadata
        metadata = metadata.split(" ") if metadata else None

    edge_list = []
    rs_meta_list = []
    if metadata:
        q_meta = getquerylist(metadata)
        rs_meta_list = filter_meta(q_meta, metadata)
    if edges:
        edge_list = getquerylist(edges)

    # get_TFlist(TFquery)
    q_tf_list, tf_name = get_tf_list(tf_query)
    rs_final = query_tf(q_tf_list, tf_name, edges, edge_list, rs_meta_list, metadata)

    # Check this:- Why there are empty spaces in my dataframe? The code works well without fillna() replacement.
    # Means there
    # are no None what I expected from mysql query. This needs to be tested.
    rs_final.replace('', np.nan, inplace=True)
    rs_final.fillna(value=np.nan, inplace=True)
    # counting the number of targets for each experiment+analysis
    # count_rs_final is a series
    count_rs_final = rs_final.count(axis=0)
    # convert series to a dict
    count_rs_finaldict = count_rs_final.to_dict()
    # convert count_rs_finaldict to a nested dict
    count_nesteddict = defaultdict(dict)
    for key, val in count_rs_finaldict.items():
        val_split = key.split('_')
        count_nesteddict['_'.join(val_split[:3])]['_'.join(val_split[3:-1])[::-1].replace('_', '.', 1)[::-1]] = val

    if not rs_final.empty:
        # if file with list of target genes is provided with -r option
        # Get the subset of results df for target genes asked in targets query file
        # Code below can handle multiple lists
        targets_mullist_dict = defaultdict(list)
        if target_genes:
            q_tg_list = []
            with open(target_genes, 'r') as q_tg:
                list_id = 'default'
                for i_q_tg in q_tg:
                    if not i_q_tg.startswith('>'):
                        q_tg_list.append(i_q_tg.strip().upper())
                        targets_mullist_dict[i_q_tg.strip().upper()].append(list_id)
                    else:
                        list_id = i_q_tg.strip('>' + whitespace)
            rs_final_res_t = rs_final[rs_final.index.isin(q_tg_list)]
        else:
            rs_final_res_t = rs_final

        # discard 'TARGET:RNASEQ'/'TARGET:CHIPSEQ' from the edges
        rs_final_res_t, chipdata_summary = tabular(rs_final_res_t)
        rs_final_trim = rs_final_res_t.replace({'TARGET:RNASEQ:': ''}, regex=True)
        rs_final_trim.replace({'INPLANTA:MICROARRAY:': ''}, regex=True, inplace=True)
        # rs_final_trim.replace('', np.nan, inplace=True)

        ############## Create a df with p-value and FC index=AGI_ID
        res_refid_dict = {}  # referenceid-metaid analysis id mapping
        refid_tf_mapping = defaultdict(list)  # referenceid-tf_name mapping
        for val_r in rs_final_trim.columns.tolist():
            res_refid_dict[int(val_r.rpartition('_')[2])] = val_r
            refid_tf_mapping[(val_r.partition('_')[0])].append('_'.join(val_r.split('_', 3)[:3]) + '_OMalleyetal_2016')

        # query the database with RefID and Ath_ID (ensures number based search)
        regulation_data = pd.DataFrame(
            Regulation.objects.select_related().filter(ref_id__in=res_refid_dict.keys()).
                values_list('ref_id', 'ath_id__agi_id', 'foldchange', 'pvalue').iterator()
            , columns=['r_refid', 'r_agiid', 'r_fc', 'r_p'])

        # filter out p-values that don't make the cutoff
        if p_value:
            regulation_data = regulation_data[regulation_data.r_p.astype(float) < p_value]

        # Note: I decided to save rs_final_trim in pickle object before adding pval and fold change.
        # As I am not going to use this in the pickle file
        # Create directory to save pickle files
        os.makedirs(output, exist_ok=True)
        rs_final_trim_output = output + '/df_jsondata.pkl'

        if not regulation_data.empty:
            regulation_data['r_p_fc'] = regulation_data['r_p'] + '||' + regulation_data['r_fc']
            regulation_data.r_refid.replace(to_replace=res_refid_dict, inplace=True)

            for name, group in regulation_data.groupby('r_refid'):
                try:
                    rs_final_trim.loc[:, name] = rs_final_trim.loc[group.r_agiid, name]
                except KeyError:
                    rs_final_trim.loc[:, name] = np.nan

            # dump rs_final_trim df. Will be used by module_createjson.
            rs_final_trim.to_pickle(rs_final_trim_output)  # dump rs_final_trim

            regulation_data_new = regulation_data.pivot(columns='r_refid', index='r_agiid', values='r_p_fc')
            # filter out rows that don't make the p-value cutoff
            # subset the df for add function. Add does not work on df with different dimensions
            regulation_data_new = regulation_data_new.loc[rs_final_trim.index]

            # code for replacing a tf name in dap_data_pivot table to experiment ids
            # if a TF has multiple experiments then it repeats the df column
            for col_trim in rs_final_trim:
                try:
                    rs_final_trim[col_trim] = rs_final_trim[col_trim].astype(str).add('||').add(
                        regulation_data_new[col_trim].values.astype(str), axis='index')
                except KeyError:
                    pass
        else:
            rs_final_trim.to_pickle(rs_final_trim_output)

        # get the dap-seq data from the database
        # dap_data = pd.DataFrame(
        #     DAPdata.objects.select_related().filter(db_tfid__db_tf_agi__in=refid_tf_mapping.keys()).values_list(
        #         'db_tfid__db_tf_agi', 'ath_id__agi_id').iterator(),
        #     columns=['DAP_tf', 'DAP_target'])
        # # dap_data.DAP_tf.replace(to_replace=refid_tf_mapping, inplace=True)
        # dap_data['present'] = 'Present'
        # dap_data_pivot = dap_data.pivot(columns='DAP_tf', index='DAP_target', values='present')

        # count_nesteddict dict will be used when creating heatmaps
        # While creating heatmaps, the dataframe is filtered based on loaded list of target genes. Means number of
        # targetgenes for
        # each exp. is a subset based on loaded target genes. Which is not a correct type2 set for hypergeometric test.
        # It should be
        # instead total number of target genes for each TF- nested dict is stored as pickles

        with open(output + '/df_eachtf_tgcount.pkl', 'wb') as pickled_totaltgs:
            pickle.dump(count_nesteddict, pickled_totaltgs, protocol=pickle.HIGHEST_PROTOCOL)

        # dump target gene lists to pickle
        target_lists = defaultdict(list)
        for key, val in targets_mullist_dict.items():
            for v in val:
                target_lists[v].append(key)

        with open(output + '/target_lists.pkl', 'wb') as f:
            pickle.dump(OrderedDict(target_lists), f, protocol=pickle.HIGHEST_PROTOCOL)

        # find an alternative for this
        # None is returned by the mysql queries and nan included by pandas
        rs_final_trim[
            rs_final_trim.apply(lambda x: x.str.match(r'^(?:None|nan)?\|\|(?:None|nan)?$').fillna(True))] = np.nan
        rs_final_trim.replace(0, np.nan, inplace=True)

        rs_final_trim.dropna(how='all', inplace=True)

        # if not dap_data_pivot.empty:
        #     dap_data_pivot_replaced = pd.concat({k: dap_data_pivot[v] for v, l in refid_tf_mapping.items() if v in
        #                                          dap_data_pivot.columns.tolist() for k in l}, axis=1)
        #     rs_reg_df_merge = rs_final_trim.merge(dap_data_pivot_replaced, how='left',
        #                                           left_index=True, right_index=True)
        # else:
        rs_reg_df_merge = rs_final_trim

        new_res_df, db_metadict, mid_tfname_dict, ath_annotation, df_genelistid = create_tabular(output,
                                                                                                 rs_reg_df_merge,
                                                                                                 chipdata_summary,
                                                                                                 targets_mullist_dict)

        # dump rs_final_res_t to a pickle object. This object will be used by create_sif to create a sif file
        rs_reg_df_merge.to_pickle(output + '/df_sif.pkl')

        out_metadata_df = getmetadata(db_metadict, mid_tfname_dict)

        # dump out_metadata_df to a pickle object. This object will be used by create_sif function to create a sif file
        out_metadata_df.to_pickle(output + '/df_metadata.pkl')

        # dump uploaded target gene list to a pickle object. This object will be used by do_clustering script
        if not df_genelistid.empty:
            df_genelistid.to_pickle(output + '/df_targetgenes.pkl')

        # dump db_tf and id_name_type dict to pickle objects. Will be used by module_createjson.
        ath_annotation.iloc[:, [1, 3]].to_pickle(output + '/df_jsonanno.pkl')

        with open(output + '/df_jsondbtf.pkl', 'wb') as pickled_json_dbtf:  # dump db_tf
            db_tf = list(TargetDBTF.objects.values_list('db_tf_agi', flat=True))  # get dbase TFs
            pickle.dump(db_tf, pickled_json_dbtf, protocol=pickle.HIGHEST_PROTOCOL)

        # return new_res_df

    else:
        raise ValueError('Query did not return any data.')

    return new_res_df, out_metadata_df


##########################################################
# Function to handle unecessary space given in the query- returns a list of elements (TF or edges)
def getquerylist(query):  # should be able to handle any user entered list: TFs or edges

    q_str = ' '.join(query).upper()

    if '[' in q_str:
        q_str = re.sub(r'[\[\]]', '', q_str.upper())

    if ' OR ' in q_str:
        q_str = q_str.replace(' OR ', ' ')

    if ' AND ' in q_str:
        q_str = q_str.replace(' AND ', ' ')

    if ' ANDNOT ' in q_str:
        q_str = q_str.replace(' ANDNOT ', ' ')

    q_list_new = list(filter(None, (x.strip() for x in q_str.split(' '))))

    return q_list_new


######################################################
# Function to convert TF query to query list: TF query could be simply a list, file or ALLTF format
def get_tf_list(tf_query):
    q_tf_list = []
    # if query is given as an input file for transcription factors
    tmptf = ' '.join(tf_query)
    # The following code can handle the queries like: -t AT2G22200 or and[tf_test.txt]
    tf_name = None
    if '.TXT' in tmptf.upper():  # if input query is a txt file
        file_index = [i for i, s in enumerate(tf_query) if '.txt' in s][0]  # get index of txt file in the list
        tf_input = tf_query[file_index].replace(']', '').replace('[', '')  # get the file name
        q_list = []
        with open(tf_input, 'r') as fl_tf:  # read the file
            for val_tf in fl_tf:
                q_list.append(val_tf.strip().upper())
        tmp_TFname = (' ' + tf_query[file_index - 1].strip().upper().replace('[', '').replace(']', '') + ' ').join(
            q_list)
        # to set the start brackets around the file elements
        s_brac = ''.join(['['] * (tf_query[file_index - 1].count('['))) + ''.join(
            ['['] * (tf_query[file_index].count('[')))
        # to set the end brackets around the file elements
        e_brac = ''.join([']'] * (tf_query[file_index].count(']')))
        my_TFname = (
            s_brac + tmp_TFname + e_brac)  # replace the file name and condition (and/or) with query constructed
        del tf_query[
            file_index - 1:file_index + 1]  # delete the file name and condition from the user provided query list
        tf_query.insert(file_index - 1,
                        my_TFname)  # insert the file name and condition with query constructed in user provided list

        tf_name = ' '.join(tf_query).split()  # split by space (bcoz the newly inserted query part is still a list)

        q_tf_list = getquerylist(tf_name)

        # print('\nFollowing is your database query:')
        # print(' '.join(tf_name))

    # if input query has all TFs: ALlTF should get all the TFs from the database
    if 'ALLTF' in tmptf.upper():
        # all_tfs = sess.query(TargetDBTF.db_tf_agi).all()
        q_tf_list = list(TargetDBTF.objects.values_list('db_tf_agi', flat=True))
        tf_name = (' ' + tf_query[0].strip().upper().replace('[', '') + ' ').join(q_tf_list).split()

    if not ('.TXT' in tmptf.upper() or 'ALLTF' in tmptf.upper()):  # if
        # input query is an expression or selection of one TF
        tf_name = [x.upper() for x in tf_query]
        q_tf_list = getquerylist(tf_name)

    return q_tf_list, tf_name


############################################################
# Function to filter pandas dataframe for user query provided
# @profile
def query_tf(q_tf_list, tf_name, edges, edgelist, rs_meta_list, metadata):
    tf_frames = []  # stores frames for all TFs
    tf_mid_map = defaultdict(list)

    for q_tf in q_tf_list:  # fetch data for each TF in a separate DF
        # print('q_tf= ',q_tf)
        edge_mid_map = defaultdict(list)
        # Combine all the TF dataframes after this loop
        tf_data = query_tfdb(q_tf, rs_meta_list)

        # if df for a Tf is empty after querying the database, don't query the DF for edges
        # or concat with other TFs data
        if not tf_data.empty:  # if tf_data df is empty, don't append to the tf_frame
            tf_edges_uniq = tf_data.EDGE.unique()
            for k, g in tf_data.groupby('EDGE'):
                edge_mid_map[k] = list(set(g['REFID']))

            for i, j in tf_data.groupby('TF'):
                refidlist_each_tf = []
                for val_test_chip in set(j['REFID']):
                    if 'CHIPSEQ' in val_test_chip:
                        refidlist_each_tf.append(val_test_chip.rpartition('_')[0])
                    else:
                        refidlist_each_tf.append(val_test_chip)
                tf_mid_map[i].extend(set(refidlist_each_tf))

            # apply and pivot_table are slower than pivot
            ## option 1
            # grouped = tf_data.groupby(['TARGET', 'REFID'], axis=0)
            # rs_gp = pd.DataFrame(grouped.EDGE.apply(lambda x: ','.join(x)).unstack('REFID'))
            ## option 2
            # rs_gp = tf_data.pivot_table(index='TARGET',columns='REFID',values='EDGE',aggfunc=lambda x: ','.join(x))
            # print('rs_gp= ',rs_gp)
            ## option 3 take only 0.1 sec compared to 7 sec with pivot_table
            # pivot does not support aggfunc
            rs_gp = tf_data.pivot(index='TARGET', columns='REFID', values='EDGE')

            if edges:  # If edges query are provided by user
                # make a TMP column in DF to tackle if an edges asked is not present in dbase for a TF
                # all rows for this column are NONE. Edge not present for a TF will look
                # for values in this column and get false as a result of expression
                rs_gp['TMP'] = None
                edgequery = create_edges_query(edges, edgelist, edge_mid_map)
                rs_gp.query(edgequery, inplace=True)  # query the df of each TF for edges
                cols_discard_rsgp = []
                for rs_gp_cols in rs_gp.columns.tolist():
                    if not rs_gp_cols in edgequery:
                        cols_discard_rsgp.append(rs_gp_cols)
                rs_gp.drop(cols_discard_rsgp, 1, inplace=True)

            rs_gp.rename(columns={x_col: x_col.rpartition('_')[0] for x_col in rs_gp.columns if 'CHIPSEQ' in x_col},
                         inplace=True)

            rs_gp_new = rs_gp.fillna('')

            if 'TMP' in rs_gp_new.columns:  # discard the tmp column from the df after query
                rs_gp_new.drop('TMP', 1, inplace=True)
            if not rs_gp_new.empty:  # if dataframe for a TF is not empty after query then append it to the multiple TFs
                tf_frames.append(rs_gp_new)  # append all the dataframes to a list

    if tf_frames:
        rs_pd_all: pd.DataFrame = pd.concat(tf_frames, axis=1,
                                            join='outer', sort=True)  # join= 'outer' represents union of multiple df
        if 'AND' in ''.join(tf_name):
            filtered_columns = rs_pd_all.columns.tolist()  # after edges were removed df contains only valid edges
            tfquery = create_tf_query(tf_name, q_tf_list, tf_mid_map, filtered_columns)
            rs_pd_all.query(tfquery, inplace=True)  # query the dataframe for intersection and complex query expression
    else:  # if no data is fetched for the given query then raise an exception and exit
        raise ValueError("No data.")
    # print('rs_pd_all.cols= ', rs_pd_all.columns.tolist())

    return rs_pd_all


################################################
# Query the database
# @profile
def query_tfdb(q_tf_name, rs_meta_list):
    rs_pd = pd.DataFrame(
        Interactions.objects.select_related().filter(db_tf_id__db_tf_agi__exact=q_tf_name).values_list(
            'db_tf_id__db_tf_agi', 'edge_id__edge_name', 'target_id__agi_id', 'ref_id__ref_id').iterator(),
        columns=['TF', 'EDGE', 'TARGET', 'REFID'])
    if rs_meta_list:
        rs_pd = rs_pd.loc[rs_pd.REFID.isin(rs_meta_list)]

    list_ref_id = rs_pd.REFID.unique()
    # print('list_ref_id= ',list_ref_id)
    meta_ref = ReferenceId.objects.select_related().filter(ref_id__in=list_ref_id).values_list('ref_id',
                                                                                               'meta_id__meta_fullid',
                                                                                               'analysis_id__analysis_fullid')

    meta_ref_dict = {val_m[0]: '{0[1]}_{0[2]}_{0[0]}'.format(val_m) for val_m in meta_ref}

    # Pandas query func throws an error if columns names are numbers so I had to include meta_id in RefID
    # column name '1', '1_2', '1_a' etc. will not work
    if not rs_pd.empty:
        rs_pd.REFID.replace(to_replace=meta_ref_dict, inplace=True)
        # pvalues '.' are replaces because pandas does not allow to use these chars with pandas.query
        rs_pd['REFID'] = rs_pd['REFID'].str.replace('.', '_')
        # attach time point with each Chipseq referenceid
        pattern = rs_pd.REFID.str.split('_').str.get(2) == 'CHIPSEQ'
        rs_pd.loc[pattern, 'REFID'] = rs_pd.loc[pattern, 'REFID'] + '_' + rs_pd.loc[pattern, 'EDGE'].str.split(
            ':').str.get(2)
    return rs_pd


##########################################################
# Function to create queries for edges
def create_edges_query(edges, edgelist, edge_mid_map):
    edgestr = ' '.join(edges)
    edges_in_mid = edgestr.upper().replace(' AND ', ' & '). \
        replace(' OR ', ' | ').replace(' ANDNOT ', ' &~ '). \
        replace('[]', 'False').replace('[', '(').replace(']', ')')
    # print 'edges_in_mid= ',edges_in_mid
    for val in edgelist:
        myval = '"%s"' % val
        count = 0
        if not len(edge_mid_map[val]) == 0:
            for i in edge_mid_map[val]:  # This code can also handle edges with multiple experiments-
                # creating OR between multiple edges- TESTING REQUIRED HERE
                if count == 0:
                    each_edge_mid = ' '.join([myval, 'in', i])
                else:
                    each_edge_mid = ' '.join(
                        [each_edge_mid, '|', myval, 'in', i])
                count += 1
            if count > 1:
                each_edge_mid = ' '.join(['(', each_edge_mid, ')'])
        if len(edge_mid_map[val]) == 0:
            each_edge_mid = ' '.join([myval, 'in', 'TMP'])
        edges_in_mid = re.sub(r'\b' + val + r'\b', each_edge_mid,
                              edges_in_mid)  # replace value in query exp with pandas query exp

    return edges_in_mid


##########################################################
# Create TF query
def create_tf_query(tf_name, q_tf_list, tf_mid_map, filtered_columns):
    tfstr = ' '.join(tf_name)
    tf_in_mid = tfstr.upper().replace(' AND ', ' & ').replace(' OR ', ' | '). \
        replace(' ANDNOT ', ' &~ ').replace('[', '(').replace(']', ')')
    for val in q_tf_list:
        count = 0
        for i in tf_mid_map[val]:
            if i in filtered_columns:
                if count == 0:
                    each_tf_mid = ' '.join([str(i) + '.notnull()'])
                else:
                    each_tf_mid = ' '.join([each_tf_mid, '|', str(i) + '.notnull()'])
            else:
                if count == 0:
                    each_tf_mid = ' '.join(['False'])
                else:
                    each_tf_mid = ' '.join([each_tf_mid, '|', 'False'])
            count += 1
        if count > 1:
            each_tf_mid = ' '.join(['(', each_tf_mid, ')'])
        tf_in_mid = tf_in_mid.replace(val, each_tf_mid)
    return tf_in_mid


def col_rename(name: str) -> str:
    if not name.endswith('_OMalleyetal_2016'):
        return "{0[0]}.{0[2]}".format(name.rpartition('_')[0].rpartition('_'))
    return name


def split_name(name: str) -> Tuple[str, str]:
    return re.match(r'^((?:[^_]*_){2}(?:[^_]*))_(.+)$', name).group(1, 2)


def combine_annotations(data, anno):
    def get_anno(col):
        try:
            return pd.concat([anno[col.name[:2]], col], ignore_index=True)
        except KeyError:
            return pd.concat([pd.Series([None] * 3), col], ignore_index=True)

    return data.apply(get_anno)


##################################
# Generate tabular output
# @profile
def create_tabular(output, rs_final_res, chipdata_summary, targets_mullist_dict):
    mid_tfname_dict = {}  # dict contains metaid to genename mapping+TF target counts
    tmp_mid_counts = {}  # dict will be used as a reference for sorting final df based on targetcounts
    mid_tfname = {}  # dict contains metaid to genename mapping
    mid_genotype_control = {}  # dict contains metaid to genotype+control mapping
    tmp_rnaseq_summary = {}
    # get the total number of genes targeted by a TF (all target genes in database)
    exp_count = querydb_exp_count(rs_final_res.columns.tolist())
    # counting number of target genes in each column (target genes as per user query)
    rs_final_res.replace('', np.nan, inplace=True)  # without this replacement it will count '' as an element
    rs_final_res.replace(0, np.nan, inplace=True)  # without this replacement it will count 0 as an element
    rs_final_res.fillna(value=np.nan, inplace=True)

    ath_annotation_query = Annotation.objects.all()

    count_series = rs_final_res.count(axis=0)
    # print('rs_final_res= ',rs_final_res)
    db_metadict = defaultdict(dict)
    # get all the metadata in a nested dict: this dict will also be used to write metadata to avoid re-query the db
    list_refid = [x.split('_')[-1] for x in rs_final_res.columns.tolist() if not x.endswith('OMalleyetal_2016')]
    # print('list_refid= ',list_refid)
    meta_info = list(Metadata.objects.select_related().filter(referenceid__ref_id__in=list_refid). \
                     values_list('meta_fullid', 'metaiddata__meta_type', 'metaiddata__meta_value'))

    for i in meta_info:
        db_metadict[i[0]][i[1]] = i[2]

    rs_final_res.dropna(how='all', inplace=True)

    ################ creating data for excel sheet headers
    for i_mid in rs_final_res.columns:
        if not i_mid.endswith('OMalleyetal_2016'):
            tf_id = i_mid.partition('_')[0]
            exp_id = '_'.join(i_mid.split('_', 3)[:3])

            tf_name = ath_annotation_query.get(agi_id=tf_id).ath_name

            # if TF does not have a name use agiid, e.g. AT2G22200- 3036 (3036)
            mid_tfname_dict[i_mid] = '{}- {} ({})'.format(tf_id if tf_name == '-' else tf_name, count_series[i_mid],
                                                          exp_count[i_mid])
            tmp_mid_counts[i_mid] = count_series[i_mid]
            mid_tfname[tf_id] = str(tf_name)

            mid_genotype_control[i_mid] = ' | '.join(db_metadict[exp_id].get(k, '') for k in
                                                     ['EXPERIMENT', 'GENOTYPE', 'TISSUE', 'CONTROL'])
            # add number of induced and repressed genes for each experiment
            if not 'CHIPSEQ' in i_mid:
                if not rs_final_res[i_mid].isnull().all():
                    induced_eachexp = (rs_final_res[i_mid].str.contains('INDUCED')).sum()
                    repressed_eachexp = (rs_final_res[i_mid].str.contains('REPRESSED')).sum()
                    tmp_rnaseq_summary[i_mid] = 'Induced-{} Repressed-{}'.format(induced_eachexp, repressed_eachexp)
                else:
                    tmp_rnaseq_summary[i_mid] = 'Induced-0 Repressed-0'

    # {**x, **y} expression is for combining two dictionaries. Here I combine chipseq and rna-seq summary
    mid_annotate_df = pd.DataFrame(
        data=[mid_tfname_dict, mid_genotype_control,
              {**chipdata_summary, **tmp_rnaseq_summary}])  # dump mid_tfname_dict to a df

    df_genelistid = pd.DataFrame.from_dict(targets_mullist_dict, orient='index')

    # Breakpoint- CODE SLOW
    df_genelistid['List___UserList'] = df_genelistid.apply(lambda x: ' '.join(filter(None, x)), axis=1)

    df_genelistid_new = df_genelistid[['List___UserList']]
    if not df_genelistid_new.empty:
        df_genelistid_new['UserList___Count'] = df_genelistid.drop('List___UserList', axis=1).count(axis=1)

    # if the target list was not given there was problem merging the empty df to the annotation df
    # if no targetgene, df=all the genes in final df
    else:
        df_genelistid_new['UserList___Count'] = np.nan
        df_genelistid_new['allgenes'] = rs_final_res.index
        df_genelistid_new.set_index(['allgenes'], inplace=True)  # set all the genes as index

    # Get the Gene names of the target genes, insert it into a df and merge with the following two dfs
    ath_annotation = pd.DataFrame(ath_annotation_query.values_list('agi_id', 'ath_name', 'ath_fullname',
                                                                   'ath_gene_type', 'ath_gene_fam').iterator(),
                                  columns=['ID___Gene ID', 'Name___Gene Name',
                                           'Full Name___Gene Full Name', 'Type___Gene Type',
                                           'Family___Gene Family']).set_index('ID___Gene ID')

    df_target_names = ath_annotation.loc[df_genelistid_new.index]
    # Remove the referenceid and replace the last occurence of '_' with '.'
    # rs_final_res.rename(columns=lambda x: x[:-2][::-1].replace('_', '.', 1)[::-1], inplace=True)

    rs_final_res.rename(columns={c: col_rename(c) for c in rs_final_res.columns}, inplace=True)
    rs_final_res.columns = pd.MultiIndex.from_tuples([split_name(c) for c in rs_final_res.columns])  # provision
    # print('Before mid_annotate_df= ',mid_annotate_df.columns.tolist())
    # Remove the referenceid and replace the last occurence of '_' with '.'
    # mid_annotate_df.rename(columns=lambda x: x[:-2][::-1].replace('_', '.', 1)[::-1], inplace=True)
    mid_annotate_df.rename(columns=col_rename, inplace=True)
    mid_annotate_df.columns = pd.MultiIndex.from_tuples([split_name(c) for c in mid_annotate_df.columns])  # provision
    # print('\n\nAfter mid_annotate_df= ',mid_annotate_df.columns.tolist())
    # print('rs_final_res= ',rs_final_res)

    # rs_final_res1 = pd.concat([mid_annotate_df, rs_final_res], axis=0)

    # print('before split res_final_res1= ',rs_final_res1.columns)

    ## multiple index code should be implemented here
    ## questionable decision
    # rs_final_res1.columns = pd.MultiIndex.from_tuples([split_name(c) for c in rs_final_res1.columns])

    # Initialize final DataFrame
    final_df = []

    for col_name, column in rs_final_res.iteritems():
        if not column.isnull().all() and column.str.contains('||', regex=False).any():
            # checks if an experiment contains fold change and pvlaues: marked with ||
            # Split 'Analysis' by || into new columns
            splitted_analysis = column.str.split('\|\|', expand=True)
            splitted_analysis.loc[(splitted_analysis[0] == ''), [1, 2]] = np.nan
            # Recreate MultiIndex
            splitted_analysis.columns = pd.MultiIndex.from_tuples(
                [(*col_name, c) for c in ['Edges', 'Pvalue', 'Log2FC']])
            # Concatenate the new columns to the final_df
            final_df.append(splitted_analysis)
        # If an experiment does not have fold change and p-value it does not split
        else:
            if col_name[1] == 'OMalleyetal_2016':
                third_level = 'DAPEdge'
            else:
                third_level = 'Edges'
            tmp_df = column.to_frame()
            tmp_df.columns = pd.MultiIndex.from_tuples([(*col_name, third_level)])
            final_df.append(tmp_df)

    final_df = pd.concat(final_df, axis=1)

    # ensure DAP column comes last
    multi_cols = []
    for name, group in groupby(final_df.columns.tolist(), key=itemgetter(0)):
        g, h = tee(group, 2)
        multi_cols.extend(filter(lambda x: x[1] != 'OMalleyetal_2016', g))
        multi_cols.extend(filter(lambda x: x[1] == 'OMalleyetal_2016', h))

    final_df = final_df[multi_cols]

    # print('final_df= ', final_df.columns)

    ############################################################

    # *****Add the gene number to the dataframe
    # rs_final_res.insert(0, 'Ind___TargetIndex', range(0, 0 + len(rs_final_res)))

    # concat with gene annotation for each target gene (by columns)
    # new_df = pd.concat([df_target_names, df_genelistid_new, rs_final_res], axis=1)
    # print('df_target_names= ',df_target_names)

    df_target_names.columns = pd.MultiIndex.from_tuples([(*split_name(c), ' ') for c in df_target_names.columns])
    df_genelistid_new.columns = pd.MultiIndex.from_tuples([(*split_name(c), ' ') for c in df_genelistid_new.columns])
    new_res_df = pd.concat([final_df, df_target_names, df_genelistid_new], axis=1, join='inner')
    # new_res_df = pd.concat([final_df, df_target_names, df_genelistid_new], axis=1)

    # concat metadata for each experiment (as headers)
    new_res_df.index.name = None  # Having index name messes with reset index
    new_res_df.reset_index(inplace=True, col_fill='GeneID')
    new_res_df.rename(columns={'index': 'ID__'}, inplace=True)
    # df_count_rows = new_res_df.shape[0]
    # new_res_df["pvalue__P"] = np.nan

    new_res_df.rename(columns=methodcaller('rstrip', '_'), inplace=True)
    new_res_df = combine_annotations(new_res_df, mid_annotate_df)

    new_res_df, total_no_exp = include_targetcount(new_res_df)  # include target count column

    # sort metaids based on number of targets hit by a TF
    mid_counts = pd.DataFrame.from_dict(tmp_mid_counts, orient='index').reset_index()
    mid_counts['index'] = mid_counts['index'].str.split('_', 3, expand=True).iloc[:, :3].apply(lambda x: '_'.join(x),
                                                                                               axis=1)
    sorted_mid_counts = mid_counts.groupby('index').sum().sort_values(0, ascending=False)

    # Change column order: Can keep the analysis of the same experiment together after sort.
    # algo is simple. sorted_mid_counts has sorted each analysis. I keep only the unique expid after sort
    # exp with max target will come first without considering the fact which analysis has more targets. It is simply
    # comparing with other experiments not within.
    multi_cols = new_res_df.columns.tolist()

    list_mid_sorted_mcols = [x_unsort for i_sort in sorted_mid_counts.index for x_unsort in multi_cols if
                             i_sort in x_unsort]

    # rearranging the columns: Columns not sorted at the moment
    multi_cols = [('Full Name', 'Gene Full Name', ' '), ('Family', 'Gene Family', ' '),
                  ('Type', 'Gene Type', ' '), ('Name', 'Gene Name', ' '), ('List', 'UserList', ' '),
                  ('UserList', 'Count', ' '), ('ID', 'GeneID', 'GeneID'), multi_cols[-1]] + list_mid_sorted_mcols
    new_res_df = new_res_df[multi_cols]

    # na_position='first' to leave the header cols (na.nan values) sorted first
    new_res_df.sort_values([('TF Count', total_no_exp)], ascending=False, inplace=True, na_position='first')

    if new_res_df.shape[0] > 1:  # Writing dataframe to excel and formatting the df excel output
        new_res_df.to_pickle(output + '/tabular_output.pkl')
    else:
        raise ValueError("Empty dataframe for your query")

    return new_res_df, db_metadict, mid_tfname_dict, ath_annotation, df_genelistid[['List___UserList']]


#########################################################################
# Counts the number of target genes in experiment and analysisIds given
def querydb_exp_count(rs_final_res_cols):
    exp_count = {}
    # print('rs_final_res_cols= ',rs_final_res_cols)
    for id_val in rs_final_res_cols:
        if not id_val.endswith('OMalleyetal_2016'):
            rs_count = Interactions.objects.filter(ref_id=id_val.split('_')[-1]).values_list(
                'target_id_id').distinct().count()
            exp_count[id_val] = int(rs_count)

    return exp_count


################################################################
# Function to include target count in dataframe
def include_targetcount(new_res_df):
    # code below is to count the Target_count: default count counts analysis
    # id for each experiment separately
    new_tmp_level_sum = new_res_df.loc[:, new_res_df.columns.get_level_values(2) == 'Edges']

    new_tmp_level_sum.replace(0, np.nan, inplace=True)
    new_tmp_level_sum.replace('', np.nan, inplace=True)

    tmp_level_sum = new_tmp_level_sum.iloc[3:, :]

    level_count = tmp_level_sum.count(axis=1)
    total_no_exp = '({})'.format(len(set(tmp_level_sum.columns.get_level_values(0))))
    # level_count contains tf count for each target (if multiple analysis count multiple time). Converting this to
    # presence/absence (1/0) at experiment id level so that tf count for an exp with multiple analysis is counted once.
    new_res_df[('TF Count', total_no_exp, '')] = level_count

    return new_res_df, total_no_exp


###################################################################
# Tabular function
def tabular(rs_final_res_t) -> Tuple[pd.DataFrame, Dict]:
    # For targetdbv2 chipseq data for one experiment (for diff time-points) is in a single column
    chipseq_cols = [col_final_df for col_final_df in rs_final_res_t if
                    'CHIPSEQ' in col_final_df]
    # In columns with chip-datasimply replacing all non digits with comma
    # and then stripping the commas at both ends of the strings
    for cs_cols in chipseq_cols:
        rs_final_res_t[cs_cols] = rs_final_res_t[cs_cols].str.replace(r"\D+", ",").str.strip(",")
    # creating summary for each ChIP-seq column
    chipdata_summary = {}
    for c_c in set(chipseq_cols):
        all_timepoints = pd.Series(
            rs_final_res_t[c_c].values.ravel()).dropna().unique().tolist()
        tmp_summary = []
        for k_c_c in list(set(all_timepoints)):
            if ',' in k_c_c:
                tmp_summary.extend([int(x_k) for x_k in k_c_c.split(',')])
            elif k_c_c:
                tmp_summary.append(int(k_c_c))
        chipdata_summary[c_c] = ':'.join(str(k) for k in sorted(set(tmp_summary)))

    for col_chip_dict in chipdata_summary:
        binary_code_dict = {}
        count_ele = len(chipdata_summary[col_chip_dict].split(':'))
        for flag, values in enumerate(chipdata_summary[col_chip_dict].split(':')):
            binary_code_dict[values] = 10 ** ((count_ele - flag) - 1)
        rs_final_res_t[col_chip_dict].replace(np.nan, '-', inplace=True)
        rs_final_res_t[col_chip_dict + '_binary'] = rs_final_res_t.apply(
            lambda x: convert_to_binary(binary_code_dict, x[col_chip_dict]),
            axis=1)
        no_zero = ''.join(['0'] * count_ele)

        rs_final_res_t[col_chip_dict] = rs_final_res_t[col_chip_dict + '_binary']

        rs_final_res_t.drop(col_chip_dict + '_binary', axis=1, inplace=True)
        rs_final_res_t[col_chip_dict].replace('-', np.nan, inplace=True)
        rs_final_res_t[col_chip_dict].replace(no_zero, np.nan, inplace=True)

    return rs_final_res_t, chipdata_summary


###################################################################
# function to convert chip-seq experiment values to binary format
def convert_to_binary(dct, entry):
    out = 0
    if not (entry == None or entry == '-'):
        for i in entry.split(','):
            if len(i) > 0:
                out += dct[i]
        return str(out).zfill(len(dct))


###################################################
# function to filter database based on metadata
def getmetadata(db_metadict, mid_tfname_dict):
    # This function takes a list of metadata ids and db_metadict as an input, write metadata to output file
    for x in mid_tfname_dict:
        db_metadict['_'.join(x.split('_')[0:3])]['*TRANSCRIPTION_FACTOR_NAME'] = mid_tfname_dict[x]
    out_metadata_df = pd.DataFrame.from_dict(data=db_metadict, orient='columns')

    return out_metadata_df


##################################
# Filter data for given METADATA
def filter_meta(q_meta, user_q_meta):
    rs_meta_tmp = []
    rs_ref = []
    rs_meta_id = []
    user_q_metastr = ' '.join(user_q_meta)
    user_q_meta_format = user_q_metastr.upper().replace(' AND ', ' & ').replace(' OR ', ' | '). \
        replace(' ANDNOT ', ' &~ ').replace('[', '(').replace(']', ')')

    # This first for loop is simply to short list the metaids that have the given value. Now it also short list the
    # metaid,
    # if an entity from the analysis is asked.
    for valm in q_meta:  # This loop is to simply get the metaids from the data for multiple conditions in the query
        # filtering the metadata and type given in the user query
        valm_split = valm.split('=')
        valm_upper = valm.upper()
        rs_meta = MetaIddata.objects.filter(meta_type__exact=valm_split[0], meta_value__exact=valm_split[1]). \
            values_list('meta_id', 'meta_type', 'meta_value')

        rs_meta_tmp.extend(x[0] for x in rs_meta)
        if 'ANALYSIS_METHOD' in valm_upper or 'ANALYSIS_CUTOFF' in valm_upper or 'ANALYSIS_BATCHEFFECT' in valm_upper:
            rs_analysis = AnalysisIddata.objects.filter(analysis_type__exact=valm_split[0],
                                                        analysis_value__exact=valm_split[1]).values_list('analysis_id',
                                                                                                         flat=True)
            rs_analysis_meta = list(ReferenceId.objects.filter(analysis_id__in=rs_analysis). \
                                    values_list('meta_id', flat=True))
            rs_meta_tmp.extend(rs_analysis_meta)

        # creating query expression- replace query with example: 'ANNA_SCHINKE in EXPERIMENTER'
        valm_format = '"{1}" in {0}'.format(*valm_split)  # create query for meta_data

        user_q_meta_format = user_q_meta_format.replace(valm, valm_format)

    # print('rs_meta_tmp= ',set(rs_meta_tmp))
    # This loop does the real job based on metaids collected from upper loop. It gets all the reference ids from meta
    #  table and
    #  all its related analysis. Combines the data from metaiddata tables and analysisiddata ito one dataframe. Now
    # pandas
    # query function is performed simply on this dataframe and list of reference ids passed the condition are
    # selected for
    # further database queries.
    # A quick note to remember- one metaid can have multiple analysis and reference id are unique for
    #  each analysis. Only a subset of reference ids will be selected for a metaid if user is making query on
    # analysis entities.
    for valm1 in set(rs_meta_tmp):  # This loop is to make combinations on the metaids identified in upper loop
        metadata_df = pd.DataFrame(list(Metadata.objects.select_related().filter(meta_id__exact=valm1). \
                                        values_list('meta_id', 'metaiddata__meta_value', 'metaiddata__meta_type',
                                                    'referenceid__ref_id',
                                                    'referenceid__analysis_id')),
                                   columns=['m_id', 'val', 'type', 'ref_id', 'analysis_id'])
        metadata_df_final = metadata_df.pivot(index='ref_id', columns='type', values='val')
        analysis_ids_meta = metadata_df.analysis_id.unique()
        analysisdata_df = pd.DataFrame()
        for vala1 in set(analysis_ids_meta):
            tmp_analysisdata_df = pd.DataFrame(list(Analysis.objects.select_related().filter(
                analysis_id__exact=vala1).values_list('analysis_id', 'analysisiddata__analysis_value',
                                                      'analysisiddata__analysis_type',
                                                      'referenceid__ref_id')),
                                               columns=['a_id', 'val', 'type', 'ref_id'])
            analysisdata_df = pd.concat([analysisdata_df, tmp_analysisdata_df])
        analysisdata_df_final = analysisdata_df.pivot(index='ref_id', columns='type', values='val')

        df_refid_query = pd.concat([metadata_df_final, analysisdata_df_final], axis=1, join='inner')

        m_df_out = df_refid_query.query(user_q_meta_format)

        if not m_df_out.empty:
            # print('final ref id passed query= ',m_df_out.index)
            rs_meta_id.extend(m_df_out.index)

    # only hits here if the query didn't return stuff
    if not rs_meta_id:
        raise ValueError('No data matched your metadata query!')

    return rs_meta_id
