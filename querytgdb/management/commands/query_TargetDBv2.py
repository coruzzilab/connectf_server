#!/usr/bin/python
# -*- coding: utf-8 -*-

from ...models import TargetDBTF, Edges, Metadata, Analysis, Annotation, ReferenceId, \
    Interactions, Regulation, MetaIddata, DAPdata
from django.core.management.base import BaseCommand, CommandError
from django.db.models import Max
from decimal import Decimal
from collections import defaultdict
import sys, os, re, operator
import numpy as np
import pandas as pd
import pickle
from querytgdb.management.commands.Modules import module_query

class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('-t', '--TFname', nargs = '+', help = 'Search by TF name or get all the data '
        'from the database (-t= OR [ALLTF]', required = True)
        parser.add_argument('-e', '--edges', nargs='+', help= 'search by edges')
        parser.add_argument('-m', '--metadata', nargs='+', help= 'search by metadata')
        parser.add_argument('-o', '--output', help= 'output file name', required=False)
        parser.add_argument('-r', '--targetgenes', help= 'list of genes provided by the user to refine the database output')

    def handle(self, *args, **options):

        self.main(options['TFname'], options['edges'], options['metadata'],
                  options['targetgenes'], options['output'])


    #######################
    # main method
    # @profile
    def main(self, TFquery, edges, metadata, targetgenes, output):

        # check if the command line arguments provided are ok
        if TFquery == None:
            raise CommandError('\nError: Either generate a table for all the TFs (--t= OR [ALLTF] \n or ' \
                               'query TFs, both can not be none \n')

        edgelist = list()
        if metadata:
            q_meta = self.getquerylist(metadata)
        if edges:
            edgelist = self.getquerylist(edges)

        self.get_TFlist(TFquery)
        q_tf_list, TFname = self.get_TFlist(TFquery)
        rs_final = self.queryTF(q_tf_list, TFname, edges, edgelist)

        if not rs_final.empty:
            # if file with list of target genes is provided with -r option
            # Get the subset of results df for target genes asked in targets query file
            # Code below can handle multiple lists
            targets_mullist_dict = defaultdict(list)
            if targetgenes:
                q_tg_list = list()
                q_tg = open(targetgenes, 'r')
                for i_q_tg in q_tg:
                    if not i_q_tg[0] == '>':
                        q_tg_list.append(i_q_tg.strip().upper())
                        targets_mullist_dict[i_q_tg.strip().upper()].append(listid)
                    else:
                        listid = i_q_tg[1:].strip()
                rs_final_res_t = rs_final[rs_final.index.isin(q_tg_list)]
            else:
                rs_final_res_t = rs_final

            # discard 'TARGET:RNASEQ'/'TARGET:CHIPSEQ' from the edges- This manipulation takes too much time (find other way)
            rs_final_trim = rs_final_res_t.groupby(rs_final_res_t.columns, axis=1). \
                apply(lambda x: x.apply(lambda y: ','.join([':'.join(l.split(':')[2:]) for l in y if pd.notnull(l)]), axis=1))
            rs_final_trim.replace('', np.nan, inplace=True)
            ############## Create a df with p-value and FC index=AGI_ID
            res_refid_dict = dict() # referenceid-metaid analysis id mapping
            refid_tf_mapping= defaultdict(list) # referenceid-TFname mapping
            for val_r in rs_final_trim.columns.tolist():
                res_refid_dict[int(val_r.split('_')[-1])] = val_r
                refid_tf_mapping[(val_r.split('_')[0])].append('_'.join(val_r.split('_')[:3])+'_OMalleyetal_2016')

            # query the database with RefID and Ath_ID (ensures number based search)
            regulation_data= pd.DataFrame(list(Regulation.objects.select_related().filter(ref_id__in=res_refid_dict.keys()).
                                          values_list('ref_id', 'ath_id__agi_id', 'foldchange', 'pvalue'))
                                          , columns=['r_refid', 'r_agiid', 'r_fc', 'r_p'])

            regulation_data['r_p_fc'] = regulation_data['r_p'] + '||' + regulation_data['r_fc']
            regulation_data.r_refid.replace(to_replace=res_refid_dict, inplace=True)
            regulation_data_new = regulation_data.pivot(columns='r_refid', index='r_agiid', values='r_p_fc')
            # subset the df for add function. Add does not work on df with different dimensions
            regulation_data_subset = regulation_data_new.loc[rs_final_trim.index]

            # get the dap-seq data from the database
            dap_data= pd.DataFrame(list(DAPdata.objects.select_related().filter(db_tfid__db_tf_agi__in=refid_tf_mapping.keys()).
                                          values_list('db_tfid__db_tf_agi', 'ath_id__agi_id'))
                                          , columns= ['DAP_tf','DAP_target'])
            dap_data.DAP_tf.replace(to_replace=refid_tf_mapping, inplace=True)
            dap_data['present']= 'Present'
            dap_data_pivot = dap_data.pivot(columns='DAP_tf',index='DAP_target',values='present')
            #print('dap_data_pivot= ',dap_data_pivot)
            dap_data= dap_data.set_index('DAP_target')
            #print('xxx index= ',dap_data.index)

            # Note: I decided to save rs_final_trim in pickle object before adding pval and fold change.
            # As I am not going to use this in the pickle file
            # Create directory to save pickle files (if the dir does not exist)
            if not os.path.exists(output + '_pickle'):  # create output directory
                os.makedirs(output + '_pickle')
            # dump rs_final_trim df. Will be used by module_createjson.
            pickled_jsondata = output + '_pickle/df_jsondata.pkl'  # dump rs_final_trim
            rs_final_trim.to_pickle(pickled_jsondata)

            #print('regulation_data_subset= ',regulation_data_subset)
            #print('rs_final_trim= ',rs_final_trim)
            #print('dap_data= ',dap_data)

            # use add func to add pvals and fold change.
            # print('rs_final_trim.shape= ',rs_final_trim.index)
            # print('regulation_data_subset.shape= ',regulation_data_subset.index)
            # ***** Make sure it adds data based on indexes (geneids) not by index position.
            rs_final_trim = rs_final_trim.astype(str).add('||')
            # print('rs_final_trim= ',rs_final_trim)
            rs_reg_df_merge_tmp = rs_final_trim.add(regulation_data_subset.values, axis='index')
            #print('rs_reg_df_merge_tmp= ',rs_reg_df_merge_tmp.shape)
            #test_df= pd.merge([rs_reg_df_merge, dap_data_pivot], how="left", left_index= True, indicator=True)
            rs_reg_df_merge = rs_reg_df_merge_tmp.merge(dap_data_pivot, how='left', left_index= True, right_index= True)
            #print('rs_reg_df_merge= ',rs_reg_df_merge)

            new_res_df, db_metadict, mid_tfname_dict, ath_annotation, df_genelistid= self.create_tabular(output, rs_reg_df_merge,
                                                                                         targetgenes, targets_mullist_dict)

            # dump rs_final_res_t to a pickle object. This object will be used by create_sif to create a sif file
            pickled_sif= output + '_pickle/df_sif.pkl'
            #rs_final_res_t.to_pickle(pickled_sif)
            rs_reg_df_merge.to_pickle(pickled_sif)

            out_metadata_df = self.getmetadata(db_metadict, mid_tfname_dict)

            # dump out_metadata_df to a pickle object. This object will be used by create_sif function to create a sif file
            pickled_metadata = output + '_pickle/df_metadata.pkl'
            out_metadata_df.to_pickle(pickled_metadata)
            # pickled_metadata.close()

            # dump uploaded target gene list to a pickle object. This object will be used by do_clustering script
            pickled_targetgenes= output + '_pickle/df_targetgenes.pkl'
            df_genelistid.to_pickle(pickled_targetgenes)

            # dump db_tf and id_name_type dict to pickle objects. Will be used by module_createjson.
            pickled_json_anno = open(output + '_pickle/df_jsonanno.pkl', 'wb')  # dump id_name_type dict
            id_name_type = {ath_x[0]: [ath_x[1], ath_x[3]] for ath_x in ath_annotation}
            pickle.dump(id_name_type, pickled_json_anno)
            pickled_json_anno.close()  # close the pickled object file

            pickled_json_dbtf = open(output + '_pickle/df_jsondbtf.pkl', 'wb')  # dump db_tf
            rs_tf = list(TargetDBTF.objects.values_list('db_tf_agi'))  # get dbase TFs
            db_tf = [x[0] for x in rs_tf]
            pickle.dump(db_tf, pickled_json_dbtf)
            pickled_json_dbtf.close()  # close the pickled object file

            return new_res_df

        else:
            message_dict = dict()
            message_dict['Warning'] = 'No Data Matched Your Query!'
            new_res_df = pd.DataFrame.from_dict(message_dict, orient='index', dtype=None)
            new_res_df.columns = ['1']
            out_metadata_df = pd.DataFrame.from_dict(message_dict, orient='index', dtype=None)
            out_metadata_df.columns = ['1']

            return new_res_df, out_metadata_df


    ##########################################################
    # Function to handle unecessary space given in the query- returns a list of elements (TF or edges)
    def getquerylist(self, query):  # should be able to handle any user entered list: TFs or edges

        q_str = ' '.join(query)

        if '[' in q_str.upper():
            q_str = q_str.upper().replace('[', '').replace(']', '')

        if ' OR ' in q_str.upper():
            q_str = q_str.upper().replace(' OR ', ' ')

        if ' AND ' in q_str.upper():
            q_str = q_str.upper().replace(' AND ', ' ')

        if ' ANDNOT ' in q_str.upper():
            q_str = q_str.upper().replace(' ANDNOT ', ' ')

        q_list_new = [x.upper().strip() for x in q_str.split(' ')]
        q_list_new = [x for x in filter(None, q_list_new)]

        return q_list_new


    ######################################################
    # Function to convert TF query to query list: TF query could be simply a list, file or ALLTF format
    def get_TFlist(self, TFquery):
        q_tf_list = list()
        # if query is given as an input file for transcription factors
        tmptf = ' '.join(TFquery)
        # The following code can handle the queries like: -t AT2G22200 or and[tf_test.txt]
        if '.TXT' in tmptf.upper():  # if input query is a txt file
            file_index = [i for i, s in enumerate(TFquery) if '.txt' in s][0]  # get index of txt file in the list
            tf_input = TFquery[file_index].replace(']', '').replace('[', '')  # get the file name
            q_list = list()
            with open(tf_input, 'r') as fl_tf:  # read the file
                for val_tf in fl_tf:
                    q_list.append(val_tf.strip().upper())
            tmp_TFname = (' ' + TFquery[file_index - 1].strip().upper().replace('[', '').replace(']', '') + ' ').join(q_list)
            # to set the start brackets around the file elements
            s_brac = ''.join(['['] * (TFquery[file_index - 1].count('['))) + ''.join(['['] * (TFquery[file_index].count('[')))
            # to set the end brackets around the file elements
            e_brac = ''.join([']'] * (TFquery[file_index].count(']')))
            my_TFname = (s_brac + tmp_TFname + e_brac)  # replace the file name and condition (and/or) with query constructed
            del TFquery[file_index - 1:file_index + 1]  # delete the file name and condition from the user provided query list
            TFquery.insert(file_index - 1,
                           my_TFname)  # insert the file name and condition with query constructed in user provided list

            TFname = ' '.join(TFquery).split()  # split by space (bcoz the newly inserted query part is still a list)

            q_tf_list = self.getquerylist(TFname)

            #print('\nFollowing is your database query:')
            #print(' '.join(TFname))

        # if input query has all TFs: ALlTF should get all the TFs from the database
        if 'ALLTF' in tmptf.upper():
            #all_tfs = sess.query(TargetDBTF.db_tf_agi).all()
            q_tf_list= list(TargetDBTF.objects.values_list('db_tf_agi', flat=True))
            TFname = (' ' + TFquery[0].strip().upper().replace('[', '') + ' ').join(q_tf_list).split()

        if not ('.TXT' in tmptf.upper() or 'ALLTF' in tmptf.upper()):  # if
            # input query is an expression or selection of one TF
            TFname = [x.upper() for x in TFquery]
            q_tf_list = self.getquerylist(TFname)

        return q_tf_list, TFname


    ############################################################
    # Function to filter pandas dataframe for user query provided
    # @profile
    def queryTF(self, q_tf_list, TFname, edges, edgelist):
        tf_frames = list()  # stores frames for all TFs
        tf_mid_map = defaultdict(list)

        for q_tf in q_tf_list:  # fetch data for each TF in a separate DF
            edge_mid_map = defaultdict(list)
            # Combine all the TF dataframes after this loop
            tf_data = module_query.queryTFDB(q_tf)

            # if df for a Tf is empty after querying the database, don't query the DF for edges
            # or concat with other TFs data
            if not tf_data.empty:  # if tf_data df is empty, don't append to the tf_frame
                tf_edges_uniq = tf_data.EDGE.unique()
                for k, g in tf_data.groupby('EDGE'):
                    edge_mid_map[k] = list(set(g['REFID']))

                for i, j in tf_data.groupby('TF'):
                    tf_mid_map[i].extend(list(set(j['REFID'])))

                # apply and pivot_table are slower than pivot
                ## option 1
                #grouped = tf_data.groupby(['TARGET', 'REFID'], axis=0)
                #rs_gp = pd.DataFrame(grouped.EDGE.apply(lambda x: ','.join(x)).unstack('REFID'))
                ## option 2
                #rs_gp = tf_data.pivot_table(index='TARGET',columns='REFID',values='EDGE',aggfunc=lambda x: ','.join(x))

                ## option 3 take only 0.1 sec compared to 7 sec with pivot_table
                rs_gp = tf_data.pivot(index='TARGET', columns='REFID', values='EDGE')

                if edges:  # If edges query are provided by user
                    # make a TMP column in DF to tackle if an edges asked is not present in dbase for a TF
                    # all rows for this column are NONE. Edge not present for a TF will look
                    # for values in this column and get false as a result of expression
                    rs_gp['TMP'] = None

                    edgequery = self.create_edges_query(edges, edgelist, edge_mid_map)
                    rs_gp.query(edgequery, inplace=True)  # query the df of each TF for edges

                if 'TMP' in rs_gp.columns:  # discard the tmp column from the df after query
                    rs_gp.drop('TMP', 1, inplace=True)
                if not rs_gp.empty:  # if dataframe for a TF is not empty after query then append it to the multiple TFs
                    tf_frames.append(rs_gp)  # append all the dataframes to a list

        if tf_frames:
            rs_pd_all = pd.concat(tf_frames, axis=1, join='outer')  # join= 'outer' represents union of multiple df
            if 'AND' in ''.join(TFname):
                filtered_columns = rs_pd_all.columns.tolist()  # after edges were removed df contains only valid edges
                tfquery = self.create_tf_query(TFname, q_tf_list, tf_mid_map, filtered_columns)
                rs_pd_all.query(tfquery, inplace=True)  # query the dataframe for intersection and complex query expression
        else:  # if no data is fetched for the given query then raise an exception and exit
            rs_pd_all = pd.DataFrame(columns=['No_data'], dtype='float')

        return rs_pd_all


    ################################################
    # Query the database
    # @profile
    def queryTFDB(self, q_TFname):

        rs= list(Interactions.objects.select_related().filter(db_tf_id__db_tf_agi__exact= q_TFname).\
                            values_list('db_tf_id__db_tf_agi','edge_id__edge_name','target_id__agi_id','ref_id__ref_id'))

        rs_pd = pd.DataFrame(rs, columns=['TF', 'EDGE', 'TARGET', 'REFID'])
        list_ref_id = rs_pd.REFID.unique()
        meta_ref= ReferenceId.objects.select_related().filter(ref_id__in= list_ref_id).\
                           values_list('ref_id','meta_id__meta_fullid','analysis_id__analysis_fullid')

        meta_ref_dict = dict()
        for val_m in meta_ref:
            meta_ref_dict[val_m[0]] = '_'.join([val_m[1], val_m[2], str(val_m[0])])

        # Pandas query func throws an error if columns names are numbers so I had to include meta_id in RefID
        # column name '1', '1_2', '1_a' etc. will not work
        if not rs_pd.empty:
            rs_pd.REFID.replace(to_replace=meta_ref_dict, inplace=True)
            # pvalues '.' are replaces because pandas does not allow to use these chars with pandas.query
            rs_pd['REFID'] = rs_pd['REFID'].str.replace('.', '_')

        return rs_pd


    ##########################################################
    # Function to create queries for edges
    def create_edges_query(self, edges, edgelist, edge_mid_map):
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
    def create_tf_query(TFname, q_tf_list, tf_mid_map, filtered_columns):
        tfstr = ' '.join(TFname)
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
                if not i in filtered_columns:
                    if count == 0:
                        each_tf_mid = ' '.join(['False'])
                    else:
                        each_tf_mid = ' '.join([each_tf_mid, '|', 'False'])
                count += 1
            if count > 1:
                each_tf_mid = ' '.join(['(', each_tf_mid, ')'])
            tf_in_mid = tf_in_mid.replace(val, each_tf_mid)
        return tf_in_mid


    ##################################
    # Generate tabular output
    #@profile
    def create_tabular(self, outfile, rs_final_res, targetgenes, targets_mullist_dict):

        mid_tfname_dict = dict()  # dict contains metaid to genename mapping+TF target counts
        tmp_mid_counts = dict()  # dict will be used as a reference for sorting final df based on targetcounts
        mid_tfname = dict()  # dict contains metaid to genename mapping
        mid_genotype_control = dict()  # dict contains metaid to genotype+control mapping
        tmp_rnaseq_summary = dict()
        # get the total number of genes targeted by a TF (all target genes in database)
        exp_count = self.querydb_exp_count(rs_final_res.columns.tolist())
        # counting number of target genes in each column (target genes as per user query)
        rs_final_res.replace('', np.nan, inplace=True)  # without this replacement it will count '' as an element
        count_series = rs_final_res.count(axis=0)
        db_metadict= defaultdict(dict)
        # get all the metadata in a nested dict: this dict will also be used to write metadata to avoid re-query the db
        list_refid= [x.split('_')[-1] for x in rs_final_res.columns.tolist() if not x.endswith('OMalleyetal_2016')]
        #print('list_refid= ',list_refid)
        meta_info= list(Metadata.objects.select_related().filter(referenceid__ref_id__in= list_refid).\
                                         values_list('meta_fullid', 'metaiddata__meta_type', 'metaiddata__meta_value'))

        for i in meta_info:
            db_metadict[i[0]][i[1]] = i[2]


        ################ creating data for excel sheet headers
        for i_mid in rs_final_res.columns:
            if not i_mid.endswith('OMalleyetal_2016'):
                tf_name= list(Annotation.objects.filter(agi_id__exact= i_mid.split('_')[0]).values_list('ath_name', flat=True))[0]

                tf_id= i_mid.split('_')[0]
                tf_exp= db_metadict['_'.join(i_mid.split('_')[0:3])]['EXPERIMENT']
                tf_control= db_metadict['_'.join(i_mid.split('_')[0:3])]['CONTROL']
                tf_geno= db_metadict['_'.join(i_mid.split('_')[0:3])]['GENOTYPE']
                if tf_name=='-': # if TF does not have a name use agiid, e.g. AT2G22200- 3036 (3036)
                    mid_tfname_dict[i_mid] = str(tf_id) + '- ' + str(count_series[i_mid]) + ' (' + str(exp_count[i_mid]) + ')'
                else:
                    mid_tfname_dict[i_mid] = str(tf_name) + '- ' + str(count_series[i_mid]) + ' (' + str(exp_count[i_mid]) + ')'
                tmp_mid_counts[i_mid] = count_series[i_mid]
                mid_tfname[tf_id] = str(tf_name)

                mid_genotype_control[i_mid] = str(tf_exp) + ' | ' + str(tf_geno) + ' | ' + str(tf_control)
                # add number of induced and repressed genes for each experiment
                if not 'CHIPSEQ' in i_mid:
                    if not rs_final_res[i_mid].isnull().all():
                        induced_eachexp = (rs_final_res[i_mid].str.contains('INDUCED') * 1).sum()
                        repressed_eachexp = (rs_final_res[i_mid].str.contains('REPRESSED') * 1).sum()
                        tmp_rnaseq_summary[i_mid] = 'Induced-' + str(induced_eachexp) + ' Repressed-' + str(repressed_eachexp)
                    else:
                        tmp_rnaseq_summary[i_mid] = 'Induced-0' + ' Repressed-0'

        # sort metaids based on number of targets hit by a TF
        #sorted_mid_counts = sorted(list(tmp_mid_counts.items()),key=operator.itemgetter(1), reverse=True)
        mid_annotate_df = pd.DataFrame(data=[mid_tfname_dict,mid_genotype_control,tmp_rnaseq_summary])  # dump mid_tfname_dict to a df
        #tmp_chip_coding = pd.DataFrame(data=chipdata_summary,index=[' '])  # dump chipdata_summary to a df

        df_genelistid = pd.DataFrame.from_dict(targets_mullist_dict, orient='index')

        # Breakpoint- CODE SLOW
        df_genelistid['List___UserList'] = df_genelistid.apply(lambda x: ' '.join(l for l in x if not l is None), axis=1)
        df_genelistid= df_genelistid['List___UserList'].to_frame()

        df_genelistid_new = df_genelistid[['List___UserList']]
        if not df_genelistid_new.empty:
            # Breakpoint- CODE SLOW
            #df_genelistid_new['UserList___Count']= df_genelistid['List___UserList'].\
            #                                        apply(lambda x: pd.value_counts(x.strip().split(' '))).sum(axis=1)
            #alternative
            df_genelistid_new['UserList___Count'] = df_genelistid['List___UserList'].str.count(' ') + 1

        # if the target list was not given there was problem merging the empty df to the annotation df
        # if no targetgene, df=all the genes in final df
        else:
            df_genelistid_new['UserList___Count'] = np.nan
            df_genelistid_new['allgenes']= rs_final_res.index
            df_genelistid_new.set_index(['allgenes'], inplace=True) # set all the genes as index

        # Get the Gene names of the target genes, insert it into a df and merge with the following two dfs
        ath_annotation= list(Annotation.objects.values_list('agi_id', 'ath_name', 'ath_fullname',
                                               'ath_gene_type', 'ath_gene_fam'))

        df_target= pd.DataFrame(ath_annotation, columns=['ID___Gene ID', 'Name___Gene Name',
                          'Full Name___Gene Full Name','Type___Gene Type','Family___Gene Family']).set_index('ID___Gene ID')
        df_target_names= df_target.loc[df_genelistid_new.index.values]
        # Remove the referenceid and replace the last occurence of '_' with '.'
        #rs_final_res.rename(columns=lambda x: x[:-2][::-1].replace('_', '.', 1)[::-1], inplace=True)


        renamedict= dict()
        for val_rename in rs_final_res.columns:
            if not val_rename.endswith('_OMalleyetal_2016'):
                renamedict[val_rename]= ('_'.join(val_rename.split('_')[:-1]))[::-1].replace('_', '.', 1)[::-1]
            else:
                renamedict[val_rename] = val_rename
        rs_final_res.rename(columns=renamedict, inplace=True)
        #print('Before mid_annotate_df= ',mid_annotate_df.columns.tolist())
        # Remove the referenceid and replace the last occurence of '_' with '.'
        #mid_annotate_df.rename(columns=lambda x: x[:-2][::-1].replace('_', '.', 1)[::-1], inplace=True)
        mid_annotate_df.rename(columns=lambda x: ('_'.join(x.split('_')[:-1]))[::-1].replace('_', '.', 1)[::-1], inplace=True)
        #print('\n\nAfter mid_annotate_df= ',mid_annotate_df.columns.tolist())
        #print('rs_final_res= ',rs_final_res)

        rs_final_res1 = pd.concat([mid_annotate_df, rs_final_res], axis=0)

        #print('before split res_final_res1= ',rs_final_res1.columns)

        ## multiple index code should be implemented here
        rs_final_res1.columns = pd.MultiIndex.from_tuples([('_'.join(c.split('_')[:3]), '_'.join(c.split('_')[3:]))
                                                        for c in rs_final_res1.columns])

        #print('after split= ',rs_final_res1.columns)

        # Initialize final DataFrame
        final_df = pd.DataFrame()

        for col_name in rs_final_res1:
            if not col_name[1]=='OMalleyetal_2016':
                if rs_final_res1[col_name].str.contains('||').all():
                    # Split 'Analysis' by || into new columns
                    splitted_analysis = rs_final_res1[col_name].str.split('\|\|', expand=True)
                    # The new column names are 0, 1, 2. Let's rename them.
                    splitted_analysis.columns = ['Edges', 'Pvalue', 'Foldchange']
                    # Recreate MultiIndex
                    splitted_analysis.columns = pd.MultiIndex.from_tuples(
                        [(col_name[0], col_name[1], c) for c in splitted_analysis.columns])
                    # Concatenate the new columns to the final_df
                    final_df = pd.concat(objs=[final_df, splitted_analysis], axis=1)
            else:
                tmp_df= pd.DataFrame(rs_final_res1[col_name])
                tmp_df.columns = pd.MultiIndex.from_tuples(
                    [(col_name[0], col_name[1], 'DAPEdge')])
                #print('tmp_df= ',tmp_df.columns)
                final_df = pd.concat(objs=[final_df, tmp_df], axis=1)

        #print('final_df= ', final_df.columns)

        ############################################################

        # *****Add the gene number to the dataframe
        #rs_final_res.insert(0, 'Ind___TargetIndex', range(0, 0 + len(rs_final_res)))

        # concat with gene annotation for each target gene (by columns)
        #new_df = pd.concat([df_target_names, df_genelistid_new, rs_final_res], axis=1)
        #print('df_target_names= ',df_target_names)

        df_target_names.columns= pd.MultiIndex.from_tuples([('_'.join(c.split('_')[:3]),'_'.join(c.split('_')[3:]), ' ')
                                                        for c in df_target_names.columns])
        df_genelistid_new.columns = pd.MultiIndex.from_tuples([('_'.join(c.split('_')[:3]), '_'.join(c.split('_')[3:]), ' ')
                                                             for c in df_genelistid_new.columns])
        new_df = pd.concat([final_df, df_target_names, df_genelistid_new], axis=1)
        #print('final_df= ',final_df.columns)
        #print('new_df= ',new_df.columns)
        #print('mid_annotate_df= ',mid_annotate_df)

        # concat metadata for each experiment (as headers)
        #new_res_df = pd.concat([mid_annotate_df, new_df], axis=0)
        new_res_df= new_df

        new_res_df.reset_index(inplace=True, col_fill='GeneID')
        new_res_df.rename(columns={'index': 'ID__'}, inplace=True)
        #df_count_rows = new_res_df.shape[0]
        #new_res_df["pvalue__P"] = np.nan

        new_res_df.rename(columns={'Full Name__': 'Full Name', 'Name__': 'Name', 'ID__': 'ID',
                        'Family__': 'Family', 'Type__': 'Type', 'List__': 'List', 'UserList__': 'UserList'}, inplace=True)

        new_res_df, total_no_exp = self.include_targetcount(new_res_df)  # include target count column

        # sort metaids based on number of targets hit by a TF
        sorted_mid_counts = sorted(list(tmp_mid_counts.items()),key=operator.itemgetter(1), reverse=True)
        # Change column order
        multi_cols = new_res_df.columns.tolist()
        #print('multi_cols= ',multi_cols)
        list_mid_aid_sorted = list(zip(*sorted_mid_counts))[0]
        #print('list_mid_aid_sorted= ',list_mid_aid_sorted)
        list_mid_sorted = [('_'.join(x.split('_')[0:3]), ('_'.join(x.split('_')[3:]))[:-2][::-1].
                            replace('_', '.', 1)[::-1]) for x in list_mid_aid_sorted]
        # rearranging the columns: Columns not sorted at the moment
        multi_cols = [('Full Name', 'Gene Full Name', ' '), ('Family', 'Gene Family', ' '),
                      ('Type', 'Gene Type', ' '), ('Name', 'Gene Name', ' '), ('List', 'UserList', ' '),
                      ('UserList', 'Count', ' '), ('ID', 'GeneID', 'GeneID')] + multi_cols[-1:] + multi_cols[1:-7]

        new_res_df = new_res_df[multi_cols]

        #na_position='first' to leave the header cols (na.nan values) sorted first
        new_res_df.sort([('Target Count', total_no_exp)], ascending=False, inplace=True, na_position='first')
        # ********************** uncomment this
        new_res_df.insert(0, ('Ind', 'Index', ''), range(-2, len(new_res_df)-2))

        #new_res_df[('Ind', 'Index','')] = np.where(new_res_df[('Ind', 'Index','')] > 0, 17)
        #new_res_df.loc[new_res_df.Ind<0] = np.nan

        # ********************** uncomment this
        #new_res_df.ix[new_res_df[('Ind', 'Index')] <=0, ('Ind', 'Index')]= np.NaN

        #print('________________________________')
        #print('new_res_df= ',new_res_df)

        if new_res_df.shape[0]> 1:  # Writing dataframe to excel and formatting the df excel output
            pk_output= outfile + '_pickle/tabular_output.pkl'
            new_res_df.to_pickle(pk_output)
        #else:
        #    print('\nNo target genes matched the query crietria!')

        return new_res_df, db_metadict, mid_tfname_dict, ath_annotation, df_genelistid


    #########################################################################
    # Counts the number of target genes in experiment and analysisIds given
    def querydb_exp_count(self, rs_final_res_cols):

        exp_count = dict()
        #print('rs_final_res_cols= ',rs_final_res_cols)
        for id_val in rs_final_res_cols:
            if not id_val.endswith('OMalleyetal_2016'):
                rs_count= Regulation.objects.filter(ref_id= id_val.split('_')[-1]).values_list('ath_id').distinct().count()
                exp_count[id_val] = int(rs_count)

        return exp_count


    ################################################################
    # Function to include target count in dataframe
    def include_targetcount(self, new_res_df):

        # code below is to count the Target_count: default count counts analysis
        # id for each experiment separately
        tmp_level_sum = (new_res_df.notnull() * 1)  # convert data to binary format to count the Target_count correctly
        tmp_level_sum.drop(['Full Name', 'Name', 'ID', 'Type', 'Family', 'List', 'UserList'],
                           axis=1, inplace=True)  # drop unecsseary columns
        tmp_level_sum.drop([0, 1, 2], axis=0, inplace=True)  # drop unecsseary rows
        tmp_level_sum.replace(0, np.nan, inplace=True)
        level_count = tmp_level_sum.sum(level=0, axis=1)
        total_no_exp = '(' + str(len(list(set(tmp_level_sum.columns.get_level_values(0))))) + ')'
        new_res_df['Target Count',total_no_exp,''] = (level_count.notnull() * 1).sum(axis=1)
        new_res_df['Target Count',total_no_exp,''] = new_res_df['Target Count',total_no_exp,''].ix[3:]. \
            astype(np.int64).astype(str)
        new_res_df['Target Count',total_no_exp,''] = new_res_df['Target Count',total_no_exp,''].ix[3:]. \
            apply(lambda x: '{: >4}'.format(x))

        return new_res_df, total_no_exp


    ###################################################
    # function to filter database based on metadata
    def getmetadata(self, db_metadict, mid_tfname_dict):
        # This function takes a list of metadata ids and db_metadict as an input, write metadata to output file
        for x in mid_tfname_dict:
            db_metadict['_'.join(x.split('_')[0:3])]['*TRANSCRIPTION_FACTOR_NAME'] = mid_tfname_dict[x]
        out_metadata_df = pd.DataFrame.from_dict(data=db_metadict, orient='columns')

        return out_metadata_df
