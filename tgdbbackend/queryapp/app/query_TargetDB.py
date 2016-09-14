#!/usr/bin/python

'''
This script can be used to query the TargetDB database. Retrieves target genes for transcription factors (TFs) in query list.
Target genes for TFs can be filtered based on edges given in query. 
This script can also be used to generate a static table for all the TFs in the TargetDB (with -e alldata option).

Input required: Database name, TF (or alldata)

Query: TFs, Edges

Output Format: Excel

Version notes: 
Current version of TargetDB database contains data for six transcription factors: 
						TGA1, NLP7, bZIP1 (WT), bZIP1 (mutant), DREB, GATA17, HSFB2A
User can construct queries for TFs, edges and metadata
All queries tested in the latest documentation runs under 20 secs
Including analysisID in this code
Last updated: Sept 6, 2016
'''

##############
# Modules
import sys, os, re
import argparse as agp
import shutil, zipfile
import numpy as np
import pandas as pd
import json
from collections import defaultdict
from sqlalchemy import create_engine,exc,func,and_,or_,distinct
from sqlalchemy.orm import sessionmaker, aliased
from create_mysqlDB import Nodes,Edges,Meta,Interactions,Genenames,Base,AccessDatabase

################################################
# Query the database
def queryTFDB(sess, q_TFname, qedgelist, qmetalist):
	
	# create alias for node1
	n1= aliased(Nodes)
	n2= aliased(Nodes) 
	count=0
	list_metaid= list()

	#print 'q_TFname= ',q_TFname
	#print 'qedgelist= ',qedgelist
	#print 'qmetalist= ',qmetalist

	if qedgelist: 
		if qmetalist: # if both edges and meta ids are given in the user query
			rs= sess.query(n1.node_name, Edges.edge_name, n2.node_name, Interactions.meta_id, Interactions.analysis_id).filter(n1.node_name==q_TFname).\
							filter(n1.node_id==Interactions.node_1_id).filter(Interactions.node_2_id==n2.node_id).\
							filter(Interactions.edge_id==Edges.edge_id).filter(Edges.edge_name.in_(qedgelist)).filter(Interactions.meta_id.in_(qmetalist)).all()
		if not qmetalist: # if only edges are given in the user query
			rs= sess.query(n1.node_name, Edges.edge_name, n2.node_name, Interactions.meta_id, Interactions.analysis_id).filter(n1.node_name==q_TFname).\
							filter(n1.node_id==Interactions.node_1_id).filter(Interactions.node_2_id==n2.node_id).\
							filter(Interactions.edge_id==Edges.edge_id).filter(Edges.edge_name.in_(qedgelist)).all()
	if qmetalist and not qedgelist: # if only metaid is given in the user query
		rs= sess.query(n1.node_name, Edges.edge_name, n2.node_name, Interactions.meta_id, Interactions.analysis_id).filter(n1.node_name==q_TFname).\
							filter(n1.node_id==Interactions.node_1_id).filter(Interactions.node_2_id==n2.node_id).\
							filter(Interactions.edge_id==Edges.edge_id).filter(Interactions.meta_id.in_(qmetalist)).all()
	if not qmetalist and not qedgelist: # if only TFs are given- no edges and no metadata asked in the query
		rs= sess.query(n1.node_name, Edges.edge_name, n2.node_name, Interactions.meta_id, Interactions.analysis_id).filter(n1.node_name==q_TFname).\
							filter(n1.node_id==Interactions.node_1_id).filter(Interactions.node_2_id==n2.node_id).\
							filter(Interactions.edge_id==Edges.edge_id).all()

	rs_pd= pd.DataFrame(rs, columns=['TF','EDGE','TARGET','META','ANALYSIS'])
#	if rs_pd.empty:
#		print 'No data matched your query!'
#		sys.exit(1)

	rs_pd['ANALYSIS'] = rs_pd['ANALYSIS'].str.replace('.','_') # pvalues '.' are replaces because pandas does not allow to use these chars with pandas.query
	rs_pd['META_ID'] = rs_pd[['META', 'ANALYSIS']].apply(lambda x: '_'.join(x), axis=1) # separating metaid and analysis id with '_' because pandas does not allow other chars
	rs_pd.drop('META', axis=1, inplace=True)
	rs_pd.drop('ANALYSIS', axis=1, inplace=True)

	if not rs_pd.empty:
		# code to change CHIPSEQ column names for time points: Example: METAID_X will be METAID_X:0, METAID_X:1, METAID_X:5
		pattern = rs_pd.META_ID.str.contains('CHIPSEQ')
		# attach time point with each metaid
		rs_pd.loc[pattern, 'META_ID'] = rs_pd.loc[pattern, 'META_ID'] + '_' + rs_pd.loc[pattern, 'EDGE'].str.split(':').str.get(2)
	return rs_pd

###################################################################################
# Function to create queries when allrnaseq, allchipseq, allinduced or allrepressed is given in the user query
def create_all_query(sess, edges, allquery, q_tf, pos):
	# This code query the database and fetches edges for a given TF (q_TF) 
	# and replaces allRNAseq (or allchipseq) with or/and between set of edges for a TF
	no1= aliased(Nodes)
	no2= aliased(Nodes)
	
	allquery_edges= list()
	tmpresult_edges= set(sess.query(Edges.edge_name).filter(no1.node_name==q_tf).\
								filter(no1.node_id==Interactions.node_1_id).filter(Interactions.node_2_id==no2.node_id).\
								filter(Interactions.edge_id==Edges.edge_id).filter(Edges.edge_name.like(allquery)).all())
	
	for x_r in tmpresult_edges:
		allquery_edges.append(x_r[0])
	all_r_edges= (' '+edges[pos-1].strip().replace('[','')+' ').join(allquery_edges)
	edges[pos]= '['+all_r_edges+']'
	return edges # return replaced edge

############################################################
# Function to filter pandas dataframe for user query provided
def queryTF(sess, q_tf_list, TFname, edges, edgelist, metalist, metadata):

	tf_frames= list() # stores frames for all TFs
	tf_mid_map= defaultdict(list)

	for q_tf in q_tf_list: # fetch data for each TF in a separate DF. Combine all the TF dataframe after this loop.
		edge_mid_map= defaultdict(list)	
		discard_pos= list()
		# This code works if ALLRNASEQ, ALLCHIPSEQ etc. is asked for edges.
		# Here I am collecting all the edges for a TF- RNAseq, chipseq data.
		# Allrnaseq or ALLCHIPSEQ will fetch edges present for each TF as opposed to looking for all dbase edges for each TF.
		all_edge_flag= 0
		if edges:
			tmp_edge= ' '.join(edges)
			if 'ALL' in tmp_edge.upper():
				all_edge_flag= all_edge_flag+1
				orig_edges= tmp_edge.split()

				for pos, val_q_ed in enumerate(edges):
					allquery= None
					if 'ALLRNASEQ' in val_q_ed.upper(): # Get ALLRNASEQ edges. This will not work if AND [ALLRNASEQ] given for edges. Go for ALLINDUCED or ALLREPRESSED in that case.
						discard_pos.append(str(pos-1))
						allquery= '%RNASEQ%'
						edges= create_all_query(sess, edges, allquery, q_tf, pos)
					if 'ALLCHIPSEQ' in val_q_ed.upper(): # Get ALLCHIPSEQ edges.
						discard_pos.append(str(pos-1))
						allquery= '%CHIPSEQ%'
						edges= create_all_query(sess, edges, allquery, q_tf, pos)
					if 'ALLINDUCED' in val_q_ed.upper(): # Get ALLINDUCED edges. This will fetch edges for pCHX as well as mCHX (if present for a TF).
						discard_pos.append(str(pos-1))
						allquery= '%INDUCED%'
						edges= create_all_query(sess, edges, allquery, q_tf, pos)
					if 'ALLREPRESSED' in val_q_ed.upper(): # Get ALLREPRESSED edges i.e. Edges for pCHX as well as mCHX.
						discard_pos.append(str(pos-1))
						allquery= '%REPRESSED%'
						edges= create_all_query(sess, edges, allquery, q_tf, pos)

				for p in sorted(discard_pos, reverse=True): # this will fail in case when both chipseq and rnaseq are provided- TEST
					del edges[int(p)] # delete all the positions with 'and' 'or' in reverse order to aviod change in the pos after deleting

				edgelist= getquerylist(edges)
		tf_data= queryTFDB(sess, q_tf, edgelist, metalist)
		
		# if df for a Tf is empty after querying the database, don't query the DF for edges or combine to the df for multiple TFs 
		if not tf_data.empty: # if tf_data df is empty, don't append to the tf_frame
			tf_edges_uniq= tf_data.EDGE.unique()

			for k,g in tf_data.groupby('EDGE'):
				edge_mid_map[k]= list(set(g['META_ID']))
			
			for i,j in tf_data.groupby('TF'): 
				tf_mid_map[i].extend(list(set(j['META_ID'])))

			grouped= tf_data.groupby(['TARGET','META_ID'], axis=0)
			rs_gp= pd.DataFrame(grouped.EDGE.apply(lambda x: ','.join(x)).unstack('META_ID'))
			
			if edges:# create query string for pandas.dataframe.query() function
				if set(edgelist)!=set(tf_edges_uniq):# Throw warning if an edge is non-existing for a TF
					print 'Warning: following edges are not present for ',q_tf,' in TargetDB:\n ',list(set(edgelist)-set(tf_edges_uniq))
					# make a TMP column in DF to tackle if an edges asked is not present in dbase for a TF
					rs_gp['TMP']= None # all rows for this column are NONE. Edge not present for a TF will look for values in this column and get false as a result of expression
				
				edgequery= create_pd_query(edges, edgelist, edge_mid_map)
				if edgequery== 'False':
					print 'No data matched your query!'
					print 'Create a new query'
					sys.exit(1)
				rs_gp.query(edgequery,inplace= True) # query the df of each TF for edges
			if 'TMP' in rs_gp.columns: # discard the tmp column from the dataframe after querying the DF
				rs_gp.drop('TMP', 1, inplace=True)
			if not rs_gp.empty: # if dataframe for a TF is not empty after quering the DF then append it to the DF for multiple TFs
				tf_frames.append(rs_gp) # append all the dataframes to a list

		if all_edge_flag!=0: # if allrnaseq or allchipseq data is asked. replace the changed edges (with actual edges name) to original edges for next TF
			edges= orig_edges
			all_edge_flag= 0

	if tf_frames:
		rs_pd_all = pd.concat(tf_frames, axis=1, join='outer') # concatenate the dfs for multiple TFs: join= 'outer' represents union of multiple dfs
	if not tf_frames: # if no data is fetched for the given query then raise an exception and exit
		print '\nNo data matched for the given query!'
		print 'Exit'
		sys.exit(1)
	
	# query df (df with multiple TFs) if intersection is asked for TFs 
	# , otherwise skip as join 'outer' used above creates union between multiple TFs
	#print '***tetsing here= ',''.join(TFname)
	if 'AND' in ''.join(TFname):
		filtered_columns= rs_pd_all.columns.tolist() # after edges were removed dataframe contains only valid edges
		tfquery= create_tf_query(TFname, q_tf_list, tf_mid_map, filtered_columns)
		#print 'tfquery= ',tfquery
		rs_pd_all.query(tfquery,inplace= True) # query the dataframe for intersection and complex query expression
	
	return rs_pd_all


##########################################################
# Create TF query
def create_tf_query(TFname, q_tf_list, tf_mid_map, filtered_columns):

	tfstr= ' '.join(TFname)
	tf_in_mid= tfstr.upper().replace(' AND ',' & ').replace(' OR ',' | ').replace(' ANDNOT ',' &~ ').replace('[','(').replace(']',')')

	for val in q_tf_list:
		count=0
		for i in tf_mid_map[val]:
			if i in filtered_columns:
				if count==0:
					each_tf_mid= ' '.join([i+'.notnull()'])
				else:
					each_tf_mid= ' '.join([each_tf_mid, '|', i+'.notnull()'])
			if not i in filtered_columns:
				if count==0:
					each_tf_mid= ' '.join(['False'])
				else:
					each_tf_mid= ' '.join([each_tf_mid, '|','False'])
			count +=1
		if count >1:
			each_tf_mid= ' '.join(['(', each_tf_mid, ')'])
		tf_in_mid= tf_in_mid.replace(val, each_tf_mid)

	return tf_in_mid


##########################################################
# Function to create queries for edges
def create_pd_query(edges, edgelist, edge_mid_map):
	
	edgestr= ' '.join(edges)
	edges_in_mid= edgestr.upper().replace(' AND ',' & ').replace(' OR ',' | ').replace(' ANDNOT ',' &~ ').\
																		replace('[]','False').replace('[','(').replace(']',')')
	#print 'edges_in_mid= ',edges_in_mid
	for val in edgelist:
		myval= '"%s"'%val
		count=0
		if not len(edge_mid_map[val])==0:
			for i in edge_mid_map[val]: # This code can also handle edges with multiple experiments- creating OR between multiple edges- TESTING REQUIRED HERE
				if count==0:
					each_edge_mid= ' '.join([myval, 'in', i])					
				else:
					each_edge_mid= ' '.join([each_edge_mid, '|', myval, 'in', i])
				count +=1
			if count >1:
				each_edge_mid= ' '.join(['(', each_edge_mid, ')'])
		if len(edge_mid_map[val])==0:
			each_edge_mid= ' '.join([myval, 'in', 'TMP'])
		edges_in_mid= re.sub(r'\b'+val+r'\b',each_edge_mid,edges_in_mid) # replace value in query expression with pandas query expressions

	return edges_in_mid


##########################################################
# Function to handle unecessary space given in the query- returns a list of elements (TF or edges)
def getquerylist(query): # should be able to handle any user entered list: TFs or edges

	q_str= ' '.join(query)
	
	if '[' in q_str.upper():
		q_str= q_str.upper().replace('[','').replace(']','')
	
	if ' OR ' in q_str.upper():
		q_str= q_str.upper().replace(' OR ',' ')
	
	if ' AND ' in q_str.upper():
		q_str= q_str.upper().replace(' AND ',' ')
		
	if ' ANDNOT ' in q_str.upper():
		q_str= q_str.upper().replace(' ANDNOT ',' ')

	q_list_new= [x.upper().strip() for x in q_str.split(' ')]
	q_list_new = filter(None, q_list_new)
	
	return q_list_new


##################################
# Filter data for given METADATA
def filter_meta(sess, q_meta, user_q_meta):

	rs_meta_tmp= list()
	rs_meta_id_tmp= list()
	user_q_metastr= ' '.join(user_q_meta)
	user_q_meta_format= user_q_metastr.upper().replace(' AND ',' & ').replace(' OR ',' | ').\
													replace(' ANDNOT ',' &~ ').replace('[','(').replace(']',')')
	for valm in q_meta: # This loop is to simply get the metaids from the data for multiple conditions in the query
		rs_meta= sess.query(Meta.meta_id).\
									filter(Meta.meta_type==valm.split('=')[0]).\
									filter(Meta.meta_value==valm.split('=')[1]) # filtering the metadata and type given in the user query
		valm_format= '"%s"'%(valm.split('=')[1])+' in '+valm.split('=')[0] # create query for meta_data
		user_q_meta_format= user_q_meta_format.replace(valm, valm_format) # creating query expression- replace query with example: 'ANNA_SCHINKE in EXPERIMENTER'
		rs_meta_tmp.extend([x[0] for x in rs_meta.all()])

	db_metadict= dict()
	for valm1 in set(rs_meta_tmp): # This loop is to make combinations on the metaids identified in upper loop
		metadata_df= pd.DataFrame(sess.query(Meta.meta_id, Meta.analysis_id, Meta.meta_value, Meta.meta_type).\
												filter(Meta.meta_id==valm1).all(), columns= ['m_id','a_id','m_val','m_type'])
		metadata_df['mid_aid']= metadata_df['m_id']+'_'+metadata_df['a_id']
		metadata_df_new= metadata_df.pivot(index='mid_aid', columns='m_type', values='m_val')
		m_df_out= metadata_df_new.query(user_q_meta_format)

		if not m_df_out.empty:			
			rs_meta_id_tmp.append(m_df_out.index[0])
	
	rs_meta_id= ['_'.join(x.split('_')[:3]) for x in rs_meta_id_tmp]

	if not rs_meta_id:
		print 'No data matched your metadata query!\n'
		sys.exit(1)

	return rs_meta_id


#########################################################################
# Counts the number of target genes in experiment and analysisIds given
def querydb_exp_count(sess, rs_final_res_cols):
	
	no1= aliased(Nodes)
	no2= aliased(Nodes)
	exp_count= dict()

	for id_val in rs_final_res_cols:
		exp_id= '_'.join(id_val.split('_')[:3]) # extract experiment id
		analys_id= '_'.join(id_val.split('_')[3:]) # extract analysis id

		rs_count= sess.query(func.count(distinct(no2.node_name))).\
				         		filter(Interactions.meta_id==exp_id).filter(Interactions.analysis_id== analys_id).\
								filter(no1.node_id==Interactions.node_1_id).filter(Interactions.node_2_id==no2.node_id).\
								filter(Interactions.edge_id==Edges.edge_id).all()
		exp_count[id_val]= int(rs_count[0][0])

	return exp_count


##################################
# Generate tabular output
def create_tabular(sess, outfile, rs_final_res, targetgenes, chipdata_summary):
	
	# rename column names- converting fdr0_01 to fdr0.01
	rs_final_res.rename(columns=lambda x: x[::-1].replace('_','.',1)[::-1], inplace=True) # replacing the last occurence of '_' with '.'
	res_df_cols= rs_final_res.columns.tolist()
	
	exp_count= querydb_exp_count(sess, res_df_cols)
	rs_final_res.replace('', np.nan, inplace=True) # without this replacement it was counting '' as an element
	count_series= rs_final_res.count(axis=0)
	# get the tf name in a dict
	mid_tfname_dict= dict() # dict will contain metaid to genename mapping+TF target counts
	mid_tfname= dict() # dict will contain metaid to genename mapping
	mid_genotype_control= dict() # dict contains metaid to gentype+control mapping
	for i_mid in rs_final_res.columns: # creating data for excel sheet headers
		tf_name= 	sess.query(Meta.meta_value).filter(Meta.meta_id==('_'.join(i_mid.split('_')[0:3]))).filter(Meta.meta_type=='TRANSCRIPTION_FACTOR_NAME').all()
		tf_id= sess.query(Meta.meta_value).filter(Meta.meta_id==('_'.join(i_mid.split('_')[0:3]))).filter(Meta.meta_type=='TRANSCRIPTION_FACTOR_ID').all()
		tf_exp= sess.query(Meta.meta_value).filter(Meta.meta_id==('_'.join(i_mid.split('_')[0:3]))).filter(Meta.meta_type=='EXPERIMENT').all()
		tf_control=	sess.query(Meta.meta_value).filter(Meta.meta_id==('_'.join(i_mid.split('_')[0:3]))).filter(Meta.meta_type=='CONTROL').all()	
		tf_geno= 	sess.query(Meta.meta_value).filter(Meta.meta_id==('_'.join(i_mid.split('_')[0:3]))).filter(Meta.meta_type=='GENOTYPE').all()
		mid_tfname_dict[i_mid]= str(tf_name[0][0])+'- '+str(count_series[i_mid])+' ('+str(exp_count[i_mid])+')'
		mid_tfname[tf_id[0][0]]= str(tf_name[0][0])
		mid_genotype_control[i_mid]= str(tf_exp[0][0])+' | '+str(tf_geno[0][0])+' | '+str(tf_control[0][0])# include the data for control here in '+' with genotype
	
	mid_tfname_df= pd.DataFrame(data=mid_tfname_dict, index=[' ']) # dump mid_tfname_dict to a df
	mid_geno_cntrl_df= pd.DataFrame(data= mid_genotype_control, index=[' ']) # dump mid_genotype_control to a df
	chip_coding= pd.DataFrame(data=chipdata_summary, index=[' ']) # dump mid_tfname_dict to a df
	chip_coding.rename(columns=lambda x: x[::-1].replace('_','.',1)[::-1], inplace=True) # replacing the last occurence of '_' with '.'
	
	# Get the Gene names of the target genes, insert it into a dataframe and merge with the following two dfs
	all_targetgenes= list(rs_final_res.index.values)
	rsall_gnames=list()

	for k in all_targetgenes:
		rs_gnames= sess.query(Genenames.ath_id, Genenames.ath_name, Genenames.ath_fullname).\
							filter(Genenames.ath_id==k.strip()).all()
		rsall_gnames.extend(rs_gnames)

	df_target_names= pd.DataFrame(rsall_gnames, columns=['ID___Gene ID','Name___Gene Name','Full Name___Gene Full Name']).set_index('ID___Gene ID')

	# concat this df to results df on metaid column names to create additional row for gene names
	new_df= pd.concat([df_target_names, rs_final_res], axis=1) # concat with gene annotation for each target gene (by columns)
	new_res_df= pd.concat([mid_tfname_df, chip_coding, new_df], axis=0) # concat metadata for each experiment (as headers)
	new_res_df.reset_index(inplace=True)

	new_res_df["pvalue___P"] = np.nan
	new_res_df.rename(columns={'index': 'ID___Gene ID'}, inplace=True)
	df_count_rows= new_res_df.shape[0]
	
	# creating multiple index by separating expid and analysisid
	new_res_df.columns= pd.MultiIndex.from_tuples([('_'.join(c.split('_')[:3]), '_'.join(c.split('_')[3:])) for c in new_res_df.columns])
	#new_res_df.rename(columns=lambda x: pd.MultiIndex.from_tuples([('_'.join(x.split('_')[:3]), '_'.join(x.split('_')[3:]))]), inplace=True) # probably a better way of renaming but not working 


	# code below is to count the Target_count: default count counts analysis id for each experiment separately
	tmp_level_sum= (new_res_df.notnull() * 1) # convert data to binary format to count the Target_count correctly
	tmp_level_sum.drop(['Full Name__','Name__','ID__','pvalue__'], axis=1, inplace=True) # drop unecsseary columns
	tmp_level_sum.drop([0,1,2], axis=0, inplace=True)# drop unecsseary rows
	tmp_level_sum.replace(0, np.nan, inplace=True)
	level_count= tmp_level_sum.sum(level=0,axis=1)
	total_no_exp= '('+str(len(list(set(tmp_level_sum.columns.get_level_values(0)))))+')'
	new_res_df['Target Count', total_no_exp]= (level_count.notnull() * 1).sum(axis=1)
	new_res_df.rename(columns={'Full Name__':'Full Name','Name__':'Name','ID__':'ID','pvalue__':'pvalue'}, inplace=True)
	multi_cols= new_res_df.columns.tolist()
	multi_cols= [('Full Name', 'Gene Full Name'), ('Name', 'Gene Name'), ('ID', 'Gene ID')]+multi_cols[-1:]+[('pvalue', 'P')]+multi_cols[1:-4] # rearraging the columns
	new_res_df= new_res_df[multi_cols]
	new_res_df.sort([('Target Count', total_no_exp)], ascending=False, inplace=True, na_position='first') # na_position='first' to leave the header cols (na.nan values) sorted first

	####################################Remove this code when implemnting 6 headers
	#new_res_df.columns = new_res_df.columns.get_level_values(0)
	new_res_df.columns = [' '.join(col).strip() for col in new_res_df.columns.values]
	new_res_df.drop(new_res_df.index[0])
	if type(new_res_df.index)==pd.MultiIndex:
		print 'new_res_df is still multiple index'
	else:
		print 'new_res_df is NOT multiple index'
	####################################

		
	
	if df_count_rows>1:# Writing dataframe to excel and formatting the excel output
		writer = pd.ExcelWriter(outfile+'/'+outfile.split('/')[-1]+'_tabular_output.xlsx') # output in excel format
		new_res_df.to_excel(writer,sheet_name='TargetDB Output') # THIS IS THE FINAL OUTPUT DATAFRAME FOR TABULAR FORMAT**
		writer= write_to_excel(writer, new_res_df)		
	else:
		print '\nNo target genes matched the query crietria!'

	return writer,new_res_df


#################################################
# function to get indegree and outdegree of a TF
# to calculat outdegree provide a transposed matrix
def get_indegree_outdegree(tf_subset_matrix_concat, TF_tofind):

	list_of_tfs= tf_subset_matrix_concat[tf_subset_matrix_concat[TF_tofind]>0].index.tolist()

	return list_of_tfs


################################################################
# function to write to excel and set the conditional formatting
def write_to_excel(writer, new_res_df):

	# get the shape of the output dataframe
	df_count_rows= new_res_df.shape[0]
	df_count_cols= new_res_df.shape[1]
	excel_count_cols= chr((df_count_cols+1) + ord('A'))

	workbook = writer.book
	worksheet = writer.sheets['TargetDB Output']
	bold_font = workbook.add_format({'bold': True, 'font_size': 13, 'border':1, 'align':'center'})
	worksheet.set_column('B:B', 15)
	worksheet.set_column('C:C', 15)
	worksheet.set_column('D:D', 15, bold_font)
	worksheet.set_column('E:E', 15, bold_font)
	worksheet.set_column('F:F', 8, bold_font)
	worksheet.set_column('G:Z', 30)
	worksheet.set_column('A:A', None, None, {'hidden': True}) # hiding meaningless column (index) created by multiindexing

	header_fmt = workbook.add_format({'font_name': 'Calibri', 'font_size': 15, 'bold': True, 'align': 'center', 'border':1})
	header_fmt_1 = workbook.add_format({'font_name': 'Calibri', 'font_size': 11, 'bold': True, 'align': 'center', 'border':1})
	worksheet.set_row(2, None, None, {'hidden': True}) # hiding unnecessary row created by multiindexing
	worksheet.set_row(3, None, header_fmt_1)
	worksheet.set_row(4, None, header_fmt)
	worksheet.set_row(5, None, header_fmt)
	format1 = workbook.add_format({'bg_color': '#FA8072', 'font_color': '#000000'})
	format2 = workbook.add_format({'bg_color': '#98FB98', 'font_color': '#000000'})
	format3 = workbook.add_format({'bg_color': '#FFFF99', 'font_color': '#000000'})
		
	# Conditonal formatting of excel sheet: Green- Induced, Red- Repressed, Yellow- CHIPSEQ
		
	worksheet.conditional_format('C6:'+excel_count_cols+str(df_count_rows+3), {'type': 'text', 'criteria': 'containing',
                                        'value': 'INDUCED', 'format': format2})
	worksheet.conditional_format('C6:'+excel_count_cols+str(df_count_rows+3), {'type': 'text', 'criteria': 'containing',
     	                                   'value': 'REPRESSED', 'format': format1})
	worksheet.conditional_format(('F7:'+excel_count_cols+str(df_count_rows+3)), {'type': 'text', 'criteria': 'containing',
                                        'value': 1,'format': format3})
	return writer


#################################
# Generate sif output
def create_sif(sess, output, tmp_df, targetgenes):

	# before creating the sif file I am combining chipseq data from different time-points (from one experiment) into one column
	# Now I am retaining only one edge for multiple time-points: e.g.: Target:CHIPSEQ:0,Target:CHIPSEQ:5 will be replaced by Target:CHIPSEQ
	sif_rs_tabular= tmp_df.groupby(tmp_df.columns, axis=1).\
							apply(lambda x:x.apply(lambda y: ','.join([l for l in y if l is not None]), axis=1))
	sif_rs_tabular.replace({'^TARGET:CHIPSEQ.*': 'TARGET:CHIPSEQ'}, regex=True, inplace=True) # replace Target:CHIPSEQ:0,Target:CHIPSEQ:5 with Target:CHIPSEQ
	sif_rs_tabular.replace('', np.nan, inplace=True)

	stacked_tmp_df= pd.DataFrame(sif_rs_tabular.stack().reset_index()) # stack converts df columns into stacked rows
	stacked_tmp_df.columns = ['TARGET','TF','EDGE'] # assign new columns
	stacked_tmp_df['TF']= stacked_tmp_df['TF'].apply(lambda x:x.split('#')[0]) # extract experimentID from experimentID+analysisID (separated by '#')
	stacked_tmp_df= stacked_tmp_df.drop_duplicates() #  dropping duplicates to solve the problem of multiple representation of same experiment (with multiple analysis ID) and chipseq multiple time-points
	stacked_tmp_df['TF']= stacked_tmp_df['TF'].apply(lambda x:x.split('_')[0]) # extract TFname from experimentID
	tf_list= list(set(stacked_tmp_df['TF'].tolist()))
	reordered_tmp_df = pd.DataFrame(stacked_tmp_df,columns=['TF','EDGE','TARGET']) # reorder the dataframe

	# SIF output in tab-delimited format
	for tf_val in tf_list:
		sub_df= reordered_tmp_df[reordered_tmp_df['TF']==tf_val]
		outfile = open(output+'/'+output.split('/')[-1]+'_'+tf_val+'.sif', 'wb') # Generates sif output file for each TF
		sub_df.to_csv(outfile,sep='\t',index=False) 
		outfile.close() # close the file resources

	outfile_all = open(output+'/'+output.split('/')[-1]+'_allTFs.sif', 'wb') # Generates sif output file for all TF
	reordered_tmp_df.to_csv(outfile_all,sep='\t',index=False)
	
	total_exp= len(sif_rs_tabular.columns.tolist()) # count the total number of experiments	
	sif_rs_tabular['target_count']= (sif_rs_tabular.notnull() * 1).sum(axis=1)
	sub_common_targets= sif_rs_tabular[sif_rs_tabular['target_count']==total_exp] # df subset of common targets (i.e. targets repersented across all the exps)
	sub_shared_targets= sif_rs_tabular[sif_rs_tabular['target_count']>1] # df subset of shared targets (i.e. target genes present in >1 exp)
	sub_common_targets.drop('target_count', 1, inplace=True) # drop the target_count column from this subset df
	sub_shared_targets.drop('target_count', 1, inplace=True) # drop the target_count column from this subset df 
	
	outfile_common = open(output+'/'+output.split('/')[-1]+'_commonTargets.sif', 'wb') # Generates sif output file for all TF
	stacked_common_targets= pd.DataFrame(sub_common_targets.stack().reset_index())
	stacked_common_targets.columns = ['TARGET','TF','EDGE']
	stacked_common_targets['TF']= stacked_common_targets['TF'].apply(lambda x:x.split('_')[0]) # extract TFname from experimentID
	stacked_common_targets[['TF','EDGE','TARGET']].to_csv(outfile_common,sep='\t',index=False)


	outfile_shared = open(output+'/'+output.split('/')[-1]+'_sharedTargets.sif', 'wb') # Generates sif output file for shared targets
	stacked_shared_targets= pd.DataFrame(sub_shared_targets.stack().reset_index())
	stacked_shared_targets.columns = ['TARGET','TF','EDGE']
	stacked_shared_targets['TF']= stacked_shared_targets['TF'].apply(lambda x:x.split('_')[0]) # extract TFname from experimentID
	stacked_shared_targets[['TF','EDGE','TARGET']].to_csv(outfile_shared,sep='\t',index=False)

	outfile_all.close() # close the file resources
	outfile_common.close() # close the file resources
	outfile_shared.close() # close the file resources

	return reordered_tmp_df	


###################################################
# function to filter database based on metadata
def getmetadata(sess, list_metaid, writer):

	# This function takes a list of metadata ids as an input, returns a nested dict that contains all the metadata types for metadata ids
	db_metadict= dict()
	expid_dict= dict()
	for valm in set(list_metaid):		
		rs_meta= sess.query(Meta.meta_id, Meta.meta_value, Meta.meta_type).filter(Meta.meta_id==valm)
		db_metadict[valm]= dict()
		for i in rs_meta:
			m_id, m_value, m_type= i
			if m_id in db_metadict:
				db_metadict[m_id][m_type]= m_value
				if m_type== 'TRANSCRIPTION_FACTOR_NAME':
					expid_dict[m_id]= m_value
	metadata_df= pd.DataFrame.from_dict(data=db_metadict,orient='columns')
	tf_name= pd.DataFrame(data=expid_dict, index=['TF NAME'])
	out_metadata_df= pd.concat([tf_name, metadata_df], axis=0) # THIS IS MY FINAL OUTPUT DATAFRAME FOR METADATA OUTPUT

	out_metadata_df.to_excel(writer,sheet_name='MetaData')
	workbook1 = writer.book
	worksheet1 = writer.sheets['MetaData']
	bold_font1 = workbook1.add_format({'bold': True, 'font_size': 13, 'border':1, 'align':'left'})
	worksheet1.set_column('A:A', 27, bold_font1)
	worksheet1.set_column('B:Z', 40)
	header_fmt = workbook1.add_format({'font_name': 'Calibri', 'font_size': 15, 'bold': True, 'align': 'center', 'border':1})
	worksheet1.set_row(1, None, header_fmt)
	writer.close()
	return out_metadata_df


############################
# main method
#@profile
def main(dbname, TFquery, edges, metadata, output, targetgenes):

	# check if the command line arguments provided are ok
	if dbname==None:
		print '\nError: Database name is not provided\n'
		sys.exit(1)
	if TFquery==None: 
		print '\nError: Either generate a table for all the TFs (--t== alltf) \n'\
			'or query based on TF (-t), both can not be none\n'
		sys.exit(1)
 
	# creating engine and session
	#engine= create_engine('mysql://coruzzilab:accesstargetdb@172.22.2.137/'+dbname)
	engine = AccessDatabase(dbname).getEngine()
	Base.metadata.bind= engine
	DBSession= sessionmaker(bind=engine)
	sess= DBSession()

	q_tf_list= list()
	rs_meta_list= list()
	edgelist= list()
	if metadata:
		q_meta= getquerylist(metadata)
		rs_meta_list= filter_meta(sess, q_meta, metadata)
	if edges:
		edgelist= getquerylist(edges)
   	
	# if query is given as an input file for transcription factors
	tmptf= ' '.join(TFquery)	
	if '.TXT' in tmptf.upper(): # if input query is a txt file
		#tf_input= TFquery[1].replace('[','').replace(']','').strip()
		file_index = [i for i, s in enumerate(TFquery) if '.txt' in s][0] # get index of txt file in the list
		tf_input= TFquery[file_index].replace(']','').replace('[','') # get the file name
		q_list= list()
		with open(tf_input, 'r') as fl_tf: # read the file
			for val_tf in fl_tf:
				q_list.append(val_tf.strip().upper())	
		tmp_TFname= (' '+TFquery[file_index-1].strip().upper().replace('[','').replace(']','')+' ').join(q_list)
		s_brac= ''.join(['[']*(TFquery[file_index-1].count('[')))+''.join(['[']*(TFquery[file_index].count('['))) # to set the start brackets around the file elements
		e_brac= ''.join([']']*(TFquery[file_index].count(']'))) # to set the end brackets around the file elements
		my_TFname= (s_brac+tmp_TFname+e_brac) # replace the file name and condition (and/or) with query constructed
		del TFquery[file_index-1:file_index+1]# delete the file name and condition from the user provided query list
		TFquery.insert(file_index-1,my_TFname)# insert the file name and condition with query constructed in user provided list
		TFname= ' '.join(TFquery).split() # split by space (bcoz the newly inserted query part is still a list)
		q_tf_list= getquerylist(TFname)
		print '\nFollowing is your database query:'
		print ' '.join(TFname)
	
	if 'ALLTF' in tmptf.upper(): # if input query has all TFs
		all_tfs= sess.query(Nodes.node_name, Nodes.node_type).filter(Nodes.node_type=='TF').all()
		for x in all_tfs:
			q_tf_list.append(x[0])
		TFname= (' '+TFquery[0].strip().upper().replace('[','')+' ').join(q_tf_list).split()

	if not ('.TXT' in tmptf.upper() or 'ALLTF' in tmptf.upper()): # if input query is an expression or selection of one TF
		TFname= [x.upper() for x in TFquery]
		q_tf_list= getquerylist(TFname)

	rs_final_res= queryTF(sess, q_tf_list, TFname, edges, edgelist, rs_meta_list, metadata)
	# if chipseq in dataframe column- exclude the time point from meta_id
	rs_final_res.columns = [str('_'.join(col.split('_')[:-1])) if 'CHIPSEQ' in col else col for col in rs_final_res.columns]
	rs_final_res.where((pd.notnull(rs_final_res)), None, inplace=True) # replacing all the Nan values in df with None

	# if file with list of target genes is provided with -r option
	# Get the subset of results dataframe (df after queries) for target genes asked in targetgenes query file
	if targetgenes:
		q_tg_list= list()
		q_tg= open(targetgenes, 'r')
		for i_q_tg in q_tg:
			q_tg_list.append(i_q_tg.strip().upper())
		rs_final_res_t= rs_final_res[rs_final_res.index.isin(q_tg_list)]
	else:
		rs_final_res_t= rs_final_res
	
	# Write Output
	if not os.path.exists(output): # create output directory
		os.makedirs(output)
	
	rs_tabular, chipdata_summary= tabular(rs_final_res_t)
	writer, new_res_df= create_tabular(sess, output, rs_tabular, targetgenes, chipdata_summary)
	reordered_tmp_df= create_sif(sess, output, rs_final_res_t, targetgenes)
	out_metadata_df= getmetadata(sess, rs_final_res.columns, writer)

	shutil.make_archive(output, 'zip', output)# create a zip file for output directory
	shutil.rmtree(output) # delete the output directory after creating zip file

	return new_res_df,out_metadata_df # returns three dfs to be displayed on user-interface


###################################################################
# Tabular function
def tabular(rs_final_res_t):

	# combine chipseq data for different time points in a single column
	rs_tabular= rs_final_res_t.groupby(rs_final_res_t.columns, axis=1).\
			apply(lambda x: x.apply(lambda y: ','.join([':'.join(l.split(':')[2:]) for l in y if l is not None]), axis=1))
	chipseq_cols= [col_final_df for col_final_df in rs_final_res_t if 'CHIPSEQ' in col_final_df]
	chipdata_summary= defaultdict(list)
	for c_c in set(chipseq_cols):
		all_timepoints= pd.Series(rs_final_res_t[c_c].values.ravel()).dropna().unique().tolist()
		chipdata_summary[c_c]= ':'.join(str(k) for k in sorted([int(x_c_dict.split(':')[2]) for x_c_dict in list(set(all_timepoints))]))

	for col_chip_dict in chipdata_summary:
		binary_code_dict= dict()
		count_ele= len(chipdata_summary[col_chip_dict].split(':'))
		for flag, values in enumerate(chipdata_summary[col_chip_dict].split(':')):					
			binary_code_dict[values]= 10**((count_ele-flag)-1)
			
		rs_tabular[col_chip_dict+'_binary'] = rs_tabular.apply(lambda x: convert_to_binary(binary_code_dict, x[col_chip_dict]), axis = 1)
		no_zero= ''.join(['0']*count_ele)
				
		rs_tabular[col_chip_dict]= 	rs_tabular[col_chip_dict+'_binary']
				
		rs_tabular.drop(col_chip_dict+'_binary', axis=1, inplace=True)
		rs_tabular[col_chip_dict].replace(no_zero,np.nan,inplace=True)			

	return rs_tabular, chipdata_summary


###################################################################
# function to convert chip-seq experiment values to binary format
def convert_to_binary(dct, entry):
	
	out = 0
	for i in entry.split(','):
		if len(i) > 0:
			out+= dct[i]  
	return str(out).zfill(len(dct))

###################################################################
# Add genesect sheet
#def get_genesect():


#	stats.hypergeom.sf(3,28774,635,137) #sf did not gave the same results as r phyper function
#	stats.hypergeom.cdf(3,28774,635,137) # got the same results as r phyper function
#	r_phyper= robjects.r['phyper']
#	r_phyper(3,635,(28774-635),137)
#	results= r_phyper(a,b,c,d)
#	results[0]
		

#############
# MAIN
#############
if __name__=='__main__':
	
	parser= agp.ArgumentParser()
	parser.add_argument('-d','--dbname',help= 'Database name',required= True)
	parser.add_argument('-t','--TFname', nargs='+', help= 'Search by TF name or'\
							 'get alldata from the database (-t alldata)', required= True)
	parser.add_argument('-e','--edges', nargs='+', help= 'Search by Edges')
	parser.add_argument('-m','--metadata', nargs='+', help= 'Search by meta data')
	parser.add_argument('-o','--output', help= 'Output file name', required= False)
	parser.add_argument('-r','--targetgenes', help= 'List of genes provided by user to refine the database output')	
	
	args= parser.parse_args()
	
	main(args.dbname, args.TFname, args.edges, args.metadata, args.output, args.targetgenes)

