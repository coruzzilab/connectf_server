'''
This module returns three json objects:
	
	database TF- database TF
	database TF- whole genome TF
	database TF- whole genome (TF+nonTF)

Notes: create_json is the master function. Returns three objects.

Last updated: September 20, 2016
'''

##############
# Modules
import sys, os
import pandas as pd
import json
import numpy as np
from create_mysqlDB import Nodes,Edges,Meta,Interactions,Genenames,Base,AccessDatabase


###################################################################
# function to create a json/object file for TFs
def create_json(sess, query_res_df, output):

	df_cols= query_res_df.columns.tolist()
	df_rows= query_res_df.index.tolist()

	rs_df_tf= sess.query(Nodes.node_name).filter(Nodes.node_type=='TF').\
							filter(Nodes.node_name.in_(df_rows)).all() # get all the TFS in the database- target genes should be in user query output dataframe
	df_tf_list= [x[0] for x in rs_df_tf] # convert database query output to list (TF list)
	rs_final_tfcolset= query_res_df[df_cols]	
		
	# function 2 to create the raw data for creating json file for whole genome TFs- Raw data means nodes, edges, size of the nodes
	output_x= output+'_genome'
	edges_json_genome, out_tf_genome, tf_genome_matrix_concat= get_json_root_edges_genome(sess, rs_final_tfcolset, df_rows)
	geneid_genome=  tf_genome_matrix_concat.columns.tolist()+tf_genome_matrix_concat.index.tolist() # For database TF matrix I need columns and rows (genome TF) as nodes
	mid_tfname_genome, tfname_tftype= id_name_mapping(sess, geneid_genome) # query the database for geneid to genename mapping
	#json_object_WG= create_json_object(geneid_genome, edges_json_genome, out_tf_genome, mid_tfname_genome, tfname_tftype, output_x)

	# function 3 to create the raw data for creating json file for TFs-Targets- Raw data means nodes, edges, size of the nodes
	output_x= output+'_targets'
	edges_json_tg, out_tg_genome, tf_tg_matrix_concat= get_json_root_edges_dbase_TGs(sess, rs_final_tfcolset, df_rows)
	geneid_TF_targets=  tf_tg_matrix_concat.columns.tolist()+tf_tg_matrix_concat.index.tolist()
	mid_tfname_dbase, tfname_tftype= id_name_mapping(sess, geneid_TF_targets) # query the database for geneid to genename mapping
	#json_object_TG= create_json_object(geneid_TF_targets, edges_json_tg, out_tg_genome, mid_tfname_dbase, tfname_tftype, output_x)
	
	# function 1 to create the raw data for creating json file for database TFs- Raw data means nodes, edges, size of the nodes
	#output_view1= output+'_dbase_view1'
	output_view1= output
	output_view2= output+'_dbase_view2'
	output_view3= output+'_dbase_view3'
	edges_json_list, out_tf_dbase, tf_subset_matrix_final= get_json_root_edges_dbase(rs_final_tfcolset, df_tf_list)
	geneid_dbase= list(set(tf_subset_matrix_final.columns.tolist()+tf_subset_matrix_final.index.tolist())) # For database TF matrix I need both rows and columns as nodes
	mid_tfname_dbase, tfname_tftype= id_name_mapping(sess, geneid_dbase) # query the database for geneid to genename mapping
	json_object_TFdbase_vw1= create_json_object(geneid_dbase, edges_json_list, out_tf_dbase, mid_tfname_dbase, tfname_tftype, output_view1)
	#json_object_TFdbase_vw2= create_json_object(geneid_dbase, edges_json_list, out_tf_genome, mid_tfname_dbase, tfname_tftype, output_view2)
	#json_object_TFdbase_vw3= create_json_object(geneid_dbase, edges_json_list, out_tg_genome, mid_tfname_dbase, tfname_tftype, output_view3)

	#return json_object_WG, json_object_TG, json_object_TFdbase_vw1, json_object_TFdbase_vw2, json_object_TFdbase_vw3
	return json_object_TFdbase_vw1

###################################################################
# Gets the root and edges for TFs present in the database
def get_json_root_edges_dbase(rs_final_tfcol_subset, df_tf_list):

	rs_final_tfsubset= rs_final_tfcol_subset[rs_final_tfcol_subset.index.isin(df_tf_list)] # get rows with TF in the database and queried by user
	tf_subset_matrix= rs_final_tfsubset.notnull() * 1 # convert edges (values in df) into binary format
	tf_subset_matrix.rename(columns=lambda x: x.split('_')[0], inplace=True) # rename the dataframe with TF name instead of experimentID
	tf_subset_matrix_concat= tf_subset_matrix.groupby(tf_subset_matrix.columns, axis=1).\
															apply(lambda p: p.apply(lambda q: sum(q.values), axis=1)) # concat dataframes for multiple TFs in dataframe columns
	uniqueele= list(set(tf_subset_matrix.columns.tolist()) - set(tf_subset_matrix.index.tolist())) # get gene id of tfs that are present in df columns but not rows
	# Replacing here all the > 0 values with 1 to count no. of genes it is each TF is targeting
	tf_subset_matrix_concat.where(tf_subset_matrix_concat==0,1,inplace=True) # Where command keeps values where condition is true and replaces false ones

	# create a dataframe for TFs missing from the rows
	if uniqueele:
		tmp_tfs_dict= dict()
		for uni_val in uniqueele:
			tmp_tfs_dict= {uni_val:[0]*len(tf_subset_matrix_concat.columns.tolist())}
		tmpdf= pd.DataFrame(tmp_tfs_dict, index=tf_subset_matrix_concat.columns.tolist()) # create a tmp df to include data for non-target dfs
		tmp_frames = [tf_subset_matrix_concat, tmpdf.transpose()]	# transpose the tmpdf as missing tfs are columns in the tmpdf df
		tf_subset_matrix_final= 	pd.concat(tmp_frames) # concat the dataframes: output+tmpdf (data for non-target dfs)
	else:
		tf_subset_matrix_final= tf_subset_matrix_concat # This is the final dataframe that contains TF-TF connections

	# unstacking the df to create edges links
	stacked_root_df = pd.DataFrame(tf_subset_matrix_final.unstack().reset_index())
	transposed_df= tf_subset_matrix_final.transpose() # transpose the df to get the outdegree- because in original df outdegree is columns
	stacked_tmp_subset_df = pd.DataFrame(transposed_df[transposed_df>0].stack().reset_index()) # make the sif like format: exclude where tfs has no connection (0 in df)
	stacked_tmp_subset_df.columns= ['source','target','val'] # assign columns to the df
	stacked_tmp_subset_df['id']= stacked_tmp_subset_df['source']+'_'+stacked_tmp_subset_df['target'] # create ID for edges (TFID_TargetID)
	stacked_tmp_subset_df.drop('val', 1, inplace=True) # we don't need the dataframe values (0 1)
	
	# merge rows of a dataframe in a dict- TF-target dict
	#print 'stacked_tmp_subset_df= ',stacked_tmp_subset_df
	edges_json_list= list() # object to be included in json final dict- that will be dupmed to a json file
	for rowdf in (stacked_tmp_subset_df.to_dict(orient='records')): # converting df to a dict to create json format list
		edges_json_dict= dict()
		edges_json_dict['data']= rowdf
		edges_json_list.append(edges_json_dict)
	
	# identify the root node: discard TF's connection to itself
	# Check point- Here the tf_subset_matrix_final will be changed- TF to itself connections will be removed to identify root node
	tf_subset_matrix_final_copy= tf_subset_matrix_final.copy(deep=True) # creating a deep copy of the df to avoid changing the df
	
	for stsd in tf_subset_matrix_final_copy:
		tf_subset_matrix_final_copy.loc[stsd, stsd] = 0 # replace values for TF-self connections to zero
	
	out_tf_dbase= (tf_subset_matrix_final_copy*1).sum(axis=0) # store this data in a dict
	tf_subset_matrix_final_copy['row_sum']= tf_subset_matrix_final_copy.sum(axis=1)
	tf_subset_matrix_final_copy['column_sum']= tf_subset_matrix_final_copy.sum(axis=0)
	#rootlist= tf_subset_matrix_final_copy[(tf_subset_matrix_final_copy['row_sum']==0) & (tf_subset_matrix_final_copy['column_sum']>0)].index.tolist()
	
	return edges_json_list, out_tf_dbase, tf_subset_matrix_final


#####################################################################################
# Gets the node sizes for TFs based on number of TFs a TF is targeting in Ath genome
def get_json_root_edges_genome(sess, rs_final_tfcol_subset, df_rows):

	rs_df_tf_genome= sess.query(Genenames.ath_id).filter(Genenames.ath_gene_type=='TXNFACTOR').\
							filter(Genenames.ath_id.in_(df_rows)).all() # get genes from df rows that are defined as TFs Ath annotation- ath_gene_type in the database
	df_tf_genome= [x[0] for x in rs_df_tf_genome] # convert database query (genome) output to list (TF list)
	rs_tfgenome_sub= rs_final_tfcol_subset[rs_final_tfcol_subset.index.isin(df_tf_genome)] # get the rows with TF (genome-wide) in database (means excluding targets that are not TFs)
	tf_genome_matrix= rs_tfgenome_sub.notnull() * 1 # convert edges (values in df) into binary format
	tf_genome_matrix.rename(columns=lambda x: x.split('_')[0], inplace=True) # rename cols (TFs+experimentID) with TFs
	tf_genome_matrix_concat= tf_genome_matrix.groupby(tf_genome_matrix.columns, axis=1).\
									apply(lambda p: p.apply(lambda q: 1 if sum(q.values) >0 else 0, axis=1)) # merge cols for tfs (1 for interaction with row TF else 0)

	# Replacing here all the > 0 values with 1 to count no. of genes it is each TF is targeting
	tf_genome_matrix_concat.where(tf_genome_matrix_concat==0,1,inplace=True) # Where keeps values where condition is true and replaces false ones
	


	# unstacking the df to create edges links
	stacked_root_df_genome = pd.DataFrame(tf_genome_matrix_concat.unstack().reset_index())
	transposed_df= tf_genome_matrix_concat.transpose() # transpose the df to get the outdegree- because in original df outdegree is columns
	stacked_tmp_subset_df_genome = pd.DataFrame(transposed_df[transposed_df>0].stack().reset_index()) # make the sif like format: exclude where tfs has no connection (0 in df)
	stacked_tmp_subset_df_genome.columns= ['source','target','val'] # assign columns to the df
	stacked_tmp_subset_df_genome['id']= stacked_tmp_subset_df_genome['source']+'_'+stacked_tmp_subset_df_genome['target'] # create ID for edges (TFID_TargetID)
	stacked_tmp_subset_df_genome.drop('val', 1, inplace=True) # we don't need the dataframe values (0 1)

	edges_json_genome= list() # object to be included in json final dict- that will be dupmed to a json file
	for rowdf in (stacked_tmp_subset_df_genome.to_dict(orient='records')): # converting df to a dict to create json format list
		edges_json_genome_dict= dict()
		edges_json_genome_dict['data']= rowdf
		edges_json_genome.append(edges_json_genome_dict)
	#print 'edges= ',edges_json_genome
	
	# Check point- It is converting all the values to binary. Be careful in case sum as data is merged from multiple cols. So far I think it is fine.
	out_tf_genome= (tf_genome_matrix_concat*1).sum(axis=0) # store this data in a dict
	#print 'out_tf_genome= ',out_tf_genome
	#print 'tf_genome_matrix_concat= ',tf_genome_matrix_concat
	return edges_json_genome, out_tf_genome, tf_genome_matrix_concat


#####################################################################################
# Gets the root and edges for TFs and Targets present in the database
def get_json_root_edges_dbase_TGs(sess, rs_final_tfcol_subset, df_tf_list):

	tf_tg_matrix= rs_final_tfcol_subset.notnull() * 1 # convert edges (values in df) into binary format
	tf_tg_matrix.rename(columns=lambda x: x.split('_')[0], inplace=True) # rename cols (TFs+experimentID) with TFs
	tf_tg_matrix_concat= tf_tg_matrix.groupby(tf_tg_matrix.columns, axis=1).\
									apply(lambda p: p.apply(lambda q: 1 if sum(q.values) >0 else 0, axis=1)) # merge cols for tfs (1 for interaction with row TF else 0)
	#print 'before= ',tf_tg_matrix_concat
	# Replacing here all the > 0 values with 1 to count no. of genes it is each TF is targeting
	tf_tg_matrix_concat.where(tf_tg_matrix_concat==0,1,inplace=True) # Where keeps values where condition is true and replaces false ones

	# unstacking the df to create edges links
	stacked_root_df_tg= pd.DataFrame(tf_tg_matrix_concat.unstack().reset_index())
	transposed_df_tg= tf_tg_matrix_concat.transpose() # transpose the df to get the outdegree- because in original df outdegree is columns
	stacked_tmp_subset_df_tg = pd.DataFrame(transposed_df_tg[transposed_df_tg>0].stack().reset_index()) # make the sif like format: exclude where tfs has no connection (0 in df)
	stacked_tmp_subset_df_tg.columns= ['source','target','val'] # assign columns to the df
	stacked_tmp_subset_df_tg['id']= stacked_tmp_subset_df_tg['source']+'_'+stacked_tmp_subset_df_tg['target'] # create ID for edges (TFID_TargetID)
	stacked_tmp_subset_df_tg.drop('val', 1, inplace=True) # we don't need the dataframe values (0 1)

	edges_json_tg= list() # object to be included in json final dict- that will be dupmed to a json file
	for rowdf in (stacked_tmp_subset_df_tg.to_dict(orient='records')): # converting df to a dict to create json format list
		edges_json_tg_dict= dict()
		edges_json_tg_dict['data']= rowdf
		edges_json_tg.append(edges_json_tg_dict)

	out_tg_genome= (tf_tg_matrix_concat*1).sum(axis=0) # store this data in a dict
	
	return edges_json_tg, out_tg_genome, tf_tg_matrix_concat


###################################################################
# create a json oject
def create_json_object(geneid_x, edges_json_list, out_tf_x, mid_tfname, tfname_tftype, output_x):

	#print 'geneid_x= ',geneid_x
	elements= dict() # CHANGE THIS TO A DICT WITH NODES AND EDGES
	# create a css dict for style sheet
	css_node= dict()
	css_node['content']= 'data(name)'
	css_node['font-family']= 'helvetica'
	css_node['font-size']= 100
	css_node['font-weight']= 'bold'
	css_node['text-valign']= 'center'
	css_node['color']= '#000000'
	css_node['shape']= 'triangle'
	css_node['background-color']= '#65CA7C'
	css_edge=dict()
	css_edge['width']= 3
	css_edge['color']= '#000000'
	css_edge['target-arrow-shape']= 'triangle'
	css_edge['background-color']= '#000000'

	nodes= list()

	for col_tfs in geneid_x:
		#print 'col_tfs= ',col_tfs
		tmp_tfdict= dict()
		tmp_tfdict['data']= dict()
		tmp_tfdict['data']['id']= col_tfs # gene ID

		if col_tfs in tfname_tftype:
			if tfname_tftype[col_tfs]=='TXNFACTOR':
				tmp_tfdict['data']['type']= tfname_tftype[col_tfs] # gene type
				tmp_tfdict['data']['color']= "#33FF52" # Assign color to TFs
				tmp_tfdict['data']['shape']= "triangle" # Assign shape to TFs
			else:
				tmp_tfdict['data']['type']= tfname_tftype[col_tfs] # gene type
				tmp_tfdict['data']['color']= "#33C1FF" # Assign color to non-TFs
				tmp_tfdict['data']['shape']= "ellipse" # Assign shape to non-TFs

		# assign node position
		tmp_tfdict['renderedPosition']= dict()
		tmp_tfdict['renderedPosition']['x']= 100
		tmp_tfdict['renderedPosition']['y']= 100
			
		# assign node size (width and height)
		tmp_tfdict['style']= dict()
		try:
			ht_wt= out_tf_x[col_tfs]
		except Exception:
			ht_wt= 0
		tmp_tfdict['style']['width']= (np.log2(ht_wt+2))*100
		tmp_tfdict['style']['height']= (np.log2(ht_wt+2))*100

		if col_tfs in mid_tfname:
			if mid_tfname[col_tfs]=='-':
				tmp_tfdict['data']['name']= col_tfs+' ('+str(int(ht_wt))+')' # gene name
			else:
				tmp_tfdict['data']['name']= mid_tfname[col_tfs]+' ('+str(int(ht_wt))+')' # gene name
		else:
			tmp_tfdict['data']['name']= col_tfs+' ('+str(int(ht_wt))+')' # gene name
			
		# assign group type
		tmp_tfdict['group']= dict()				
		tmp_tfdict['group']= 'nodes'
			
		nodes.append(tmp_tfdict) # append the style dict to a element list
	
	elements['nodes']= nodes	
	elements['edges']= edges_json_list

	json_output_dict= dict()
	json_output_dict['elements']= dict() # create elements dict
	json_output_dict['elements']= elements
	json_output_dict['layout']= dict() # Set layout
	json_output_dict['layout']['name']= 'preset'
	json_object= json.dumps(json_output_dict)
 
	if output_x:
		with open(output_x.split('/')[-1]+'_cy.json', 'wb') as out_jsonfile:
			json.dump(json_output_dict, out_jsonfile, sort_keys = True, indent = 4, ensure_ascii=False)

	return json_object


###################################################################
# Function to create gene_id to gene_name mapping
def id_name_mapping(sess, geneidlist):

	mid_tfname= dict()
	tfname_tftype= dict()
	for i_mid in geneidlist:

		tf_name= 	sess.query(Genenames.ath_id, Genenames.ath_name, Genenames.ath_gene_type).\
									filter(Genenames.ath_id==i_mid).all()
		if 	tf_name:						
			mid_tfname[tf_name[0][0]]= str(tf_name[0][1])
			tfname_tftype[tf_name[0][0]]= str(tf_name[0][2])
		else:
			mid_tfname[i_mid]= '-'
			tfname_tftype[i_mid]= 'unknown'

	return mid_tfname, tfname_tftype


