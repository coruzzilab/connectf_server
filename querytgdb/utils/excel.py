import os
import shutil
from collections import defaultdict

import pandas as pd

__all__ = ('create_excel_zip',)


def create_excel_zip(pickledir):
    writer, outdirpath = read_pickled_targetdbout(pickledir)  # Func. to write
    read_pickled_metadata(pickledir, writer) # writer allows it to write to the same excel different sheet
    return create_sif(pickledir, outdirpath)


########################################################################################
# function to read TargetDB output dataframe from pickles and writes to an excel sheet
def read_pickled_targetdbout(pickledir):
    pickled_pandas = pd.read_pickle(pickledir + '/' + 'tabular_output.pkl')

    # create an output directory for downloadable zip file
    outdir = pickledir.replace('_pickle', '')
    if not os.path.exists(outdir):  # create output directory
        os.makedirs(outdir)
    else:
        shutil.rmtree(outdir)
        os.makedirs(outdir)

    # get the absolute path of the directory
    outdirpath = os.path.abspath(outdir)

    writer = pd.ExcelWriter(outdir + '/tabular_output.xlsx',
                            engine='xlsxwriter')  # output in excel format
    pickled_pandas.to_excel(writer, sheet_name='TargetDB Output')  # FINAL OUTPUT DATAFRAME

    # FOR TABULAR FORMAT**
    writer = write_to_excel(writer, pickled_pandas)

    return writer, outdirpath


################################################################
# function to write to excel and set the conditional formatting
def write_to_excel(writer, new_res_df):
    # get the shape of the output dataframe
    df_count_rows = new_res_df.shape[0]
    df_count_cols = new_res_df.shape[1]

    header_dict = defaultdict(list)
    #id_headerdict = defaultdict(list)
    id_headerdict= dict()
    header_df = new_res_df.iloc[0:3]
    header_df.drop(['Full Name', 'Name', 'ID', 'Type', 'Family', 'List', 'UserList', 'TF Count'],
                   axis=1, inplace=True)  # drop unecsseary columns

    #####################################################################################
    # code to get which columns needs to be merged for the header annotation
    ## Can be moved to a different function
    listlen = 0
    header_list = header_df.columns.tolist()
    count = 10 # start merging with column 10/J. Columns before J are annotation columns
    for val_hls_tmp in header_list:
        id_headerdict[val_hls_tmp[0]] = defaultdict(list)
    for val_hls in header_list:
        listlen = listlen + 1
        lines = None
        if count == 10:
            prev_count = 10
            prev_header_level0 = val_hls[0]
            prev_header_level1 = val_hls[1]
        if val_hls[2] == 'Edges':
            id_headerdict[val_hls[0]][val_hls[1]].append(header_df.xs(0, drop_level=True)[(val_hls[0], val_hls[1], 'Edges')])
            id_headerdict[val_hls[0]][val_hls[1]].append(header_df.xs(1, drop_level=True)[(val_hls[0], val_hls[1], 'Edges')])
            id_headerdict[val_hls[0]][val_hls[1]].append(header_df.xs(2, drop_level=True)[(val_hls[0], val_hls[1], 'Edges')])
        if val_hls[0] == prev_header_level0 and val_hls[1] == prev_header_level1:
            count += 1
            if listlen == len(header_list):
                lines = str(prev_count) + ':' + str(count - 1)
                header_dict[lines].append(id_headerdict[prev_header_level0][prev_header_level1][0])  # get annotation at index 0
                header_dict[lines].append(id_headerdict[prev_header_level0][prev_header_level1][1])  # get annotation at index 1
                header_dict[lines].append(id_headerdict[prev_header_level0][prev_header_level1][2])  # get annotation at index 2
        else:
            if not prev_header_level1 == 'OMalleyetal_2016':
                lines = str(prev_count) + ':' + str(count - 1)
                if val_hls[1] == 'OMalleyetal_2016':
                    count = count + 1
                    prev_count = count
                else:
                    prev_count = count
                    count = count + 1
                if lines:
                    header_dict[lines].append(id_headerdict[prev_header_level0][prev_header_level1][0])  # get annotation at index 0
                    header_dict[lines].append(id_headerdict[prev_header_level0][prev_header_level1][1])  # get annotation at index 1
                    header_dict[lines].append(id_headerdict[prev_header_level0][prev_header_level1][2])  # get annotation at index 2
                prev_header_level0 = val_hls[0]
                prev_header_level1 = val_hls[1]
            else:
                prev_header_level0 = val_hls[0]
                prev_header_level1 = val_hls[1]
                prev_count = count
                count = count + 1
    #####################################################################################

    excel_count_cols = colToExcel(df_count_cols + 1)

    workbook = writer.book
    worksheet = writer.sheets['TargetDB Output']
    bold_font = workbook.add_format({'bold': True, 'font_size': 13, 'border': 1, 'align': 'center'})
    align_font = workbook.add_format({'font_size': 13, 'align': 'center'})
    worksheet.set_column('B:F', 15)
    worksheet.set_column('G:G', 10, align_font)
    worksheet.set_column('H:I', 15, bold_font)
    worksheet.set_column('J:' + excel_count_cols, 15)
    worksheet.set_column('A:A', None, None, {'hidden': True})  # hiding meaningless column created by multiindexing

    header_fmt = workbook.add_format({'font_name': 'Calibri', 'font_size': 14,
                                      'bold': True, 'align': 'center', 'border': 1})
    header_fmt_1 = workbook.add_format({'font_name': 'Calibri', 'font_size': 11,
                                        'bold': True, 'align': 'center', 'border': 1})
    # worksheet.set_row(2, None, None, {'hidden': True})  # hiding unnecessary row created by multiindexing
    worksheet.set_row(3, None, header_fmt)
    worksheet.set_row(4, None, header_fmt)
    worksheet.set_row(5, None, header_fmt)
    worksheet.set_row(6, None, header_fmt)
    format1 = workbook.add_format({'bg_color': '#FA8072', 'font_color': '#000000'})
    format2 = workbook.add_format({'bg_color': '#98FB98', 'font_color': '#000000'})
    format3 = workbook.add_format({'bg_color': '#FFFF99', 'font_color': '#000000'})
    format4 = workbook.add_format({'bg_color': '#F7DC6F', 'font_color': '#000000'})

    # Merge format
    merge_format = workbook.add_format({
        'bold': 1,
        'font_size': 14,
        'border': 1,
        'align': 'center',
        'valign': 'vcenter',
        'fg_color': 'white'})

    # Merging columns
    for val_merge in header_dict:
        merge_start = colToExcel(int(val_merge.split(':')[0]))
        merge_end = colToExcel(int(val_merge.split(':')[1]))
        # print('merge_start= ',(merge_start+'5:'+merge_end+'5'))
        worksheet.merge_range(merge_start + '4:' + merge_end + '4', '', merge_format)
        worksheet.merge_range(merge_start + '5:' + merge_end + '5', header_dict[val_merge][0], merge_format)
        worksheet.merge_range(merge_start + '6:' + merge_end + '6', header_dict[val_merge][1], merge_format)
        worksheet.merge_range(merge_start + '7:' + merge_end + '7', header_dict[val_merge][2], merge_format)

    # Conditonal formatting of excel sheet: Green- Induced, Red- Repressed, Yellow- CHIPSEQ
    worksheet.conditional_format('J8:' + excel_count_cols + str(df_count_rows + 3),
                                 {'type': 'text', 'criteria': 'containing', 'value': 'INDUCED', 'format': format2})
    worksheet.conditional_format('J8:' + excel_count_cols + str(df_count_rows + 3),
                                 {'type': 'text', 'criteria': 'containing', 'value': 'REPRESSED', 'format': format1})
    # worksheet.conditional_format(
    #    ('K7:' + excel_count_cols + str(df_count_rows + 3)),
    #    {'type': 'text', 'criteria': 'containing',
    #     'value': 1, 'format': format3})
    worksheet.conditional_format('J8:' + excel_count_cols + str(df_count_rows + 3),
                                 {'type': 'text', 'criteria': 'containing', 'value': 'Present', 'format': format4})

    return writer


#################################
# get excel style column name
def colToExcel(col):  # col is 1 based
    excelCol = str()
    div = col
    while div:
        (div, mod) = divmod(div - 1, 26)  # will return (x, 0 .. 25)
        excelCol = chr(mod + 65) + excelCol

    return excelCol


################################
# Write metadata to excel
def read_pickled_metadata(pickledir, writer):
    # Read Metadata pickle
    pickled_metadf = pd.read_pickle(pickledir + '/' + 'df_metadata.pkl')

    # Write to excel
    pickled_metadf.to_excel(writer, sheet_name='MetaData')
    workbook1 = writer.book
    worksheet1 = writer.sheets['MetaData']
    bold_font1 = workbook1.add_format({'bold': True, 'font_size': 13, 'border': 1, 'align': 'left'})
    worksheet1.set_column('A:A', 27, bold_font1)
    excel_count_cols = colToExcel((pickled_metadf.shape[1]) + 1)
    worksheet1.set_column('B:' + excel_count_cols, 40)
    header_fmt = workbook1.add_format(
        {'font_name': 'Calibri', 'font_size': 15, 'bold': True, 'align': 'center', 'border': 1})
    worksheet1.set_row(1, None, header_fmt)
    writer.close()


#################################
# Generate sif output
# @profile
def create_sif(pickledir, output):
    # read the pickled dataframe
    pickled_sifdf_tmp = pd.read_pickle(pickledir + '/' + 'df_sif.pkl')

    # Before creating the sif file exclude the dap-seq column. If required later, find a way to add the dap-seq
    # validation column to the .tbl files
    keep_cols = list()
    for val_col in pickled_sifdf_tmp.columns:
        if not val_col.endswith('_OMalleyetal_2016'):
            keep_cols.append(val_col)

    pickled_sifdf = pickled_sifdf_tmp[keep_cols]

    # ** Very Slow (10 sec) Not sure if I need this anymore
    sif_rs_tabular = pickled_sifdf.groupby(pickled_sifdf.columns, axis=1). \
        apply(lambda x: x.apply(lambda y: ','.join([l for l in y if pd.notnull(l)]), axis=1))

    sif_rs_tabular.replace('', pd.np.nan, inplace=True)  # purpose of replacing spaces with nan??

    stacked_tmp_df = pd.DataFrame(
        sif_rs_tabular.stack().reset_index())  # stack converts df columns into stacked rows
    stacked_tmp_df.columns = ['TARGET', 'TF', 'EDGE']  # assign new columns

    # extract experimentID from experimentID_analysisID
    stacked_tmp_df['TF'] = stacked_tmp_df['TF'].apply(lambda x: '_'.join(x.split('_')[:3]))
    stacked_tmp_df = stacked_tmp_df.drop_duplicates()  # dropping duplicates to solve the problem of multiple
    # representation of same experiment (with multiple analysis ID) and chipseq multiple time-points
    stacked_tmp_df['TF'] = stacked_tmp_df['TF'].apply(lambda x: x.split('_')[0])  # extract TFname from expID

    tf_list = list(set(stacked_tmp_df['TF'].tolist()))
    reordered_tmp_df = pd.DataFrame(stacked_tmp_df, columns=['TF', 'EDGE', 'TARGET'])  # reorder the df

    # following is separating pvalue and fold-change from EDGE column.
    # Regular expression parses the code and assigns the column names at the same time
    # For example '(?P<Col2>.*)\|{2,}' will grab everything up to the first double | and call it Col2
    # regex = '(?P<EDGE>.*)\|{2,}(?P<PVALUE>.*)\((?P<FOLDCHANGE>.*)\)'

    # First for chipseq or any other data where edge does not have a pval and fc. Simply add the pattern '||-||-'
    # this will reflect pval and fc with -
    reordered_tmp_df.EDGE= reordered_tmp_df.EDGE.apply(lambda x: x+'||-||-' if not '||' in x else x)
    regex = '(?P<EDGE>.*)\|{2,}(?P<PVALUE>.*)\|{2,}(?P<FOLDCHANGE>.*)'
    reordered_tmp_df = reordered_tmp_df.assign(**reordered_tmp_df.EDGE.str.
                                               extract(regex, expand=True).to_dict('list'))

    # regular expression '^[01]+$' denotes that the string should only contain binary number
    # ^ means start with 0 or 1, [01]+ contains only 0 and 1, $ means ends with 0 or 1
    reordered_tmp_df['EDGE'].replace(to_replace='^[01]+$', value='CHIPSEQ', regex=True, inplace=True)

    #reordered_tmp_df = reordered_tmp_df.loc[reordered_tmp_df.EDGE != '']  # remove the edges that are not linked
    # ALTERNATIVE
    # reordered_tmp_df.EDGE.str.extract(regex, expand=True).combine_first(reordered_tmp_df)

    outfile_all = open(output + '/' + output.split('/')[-1] + '_allTFs.sif',
                       'w')  # Generates sif output file for all TF
    outfile_all_tbl = open(output + '/' + output.split('/')[-1] + '_allTFs.tbl',
                           'w')  # Generates sif output file for all TF
    reordered_tmp_df[['TF', 'EDGE', 'TARGET']].to_csv(outfile_all, sep='\t', index=False, header=False)
    reordered_tmp_df['shared name'] = reordered_tmp_df['TF'] + ' (' + reordered_tmp_df['EDGE'] + ') ' + \
                                      reordered_tmp_df['TARGET']
    reordered_tmp_df[['shared name', 'FOLDCHANGE', 'PVALUE']].to_csv(outfile_all_tbl, sep='\t', index=False)

    create_genelists(reordered_tmp_df, output)

    # SIF output in tab-delimited format
    for tf_val in tf_list:
        sub_df = reordered_tmp_df[reordered_tmp_df['TF'] == tf_val]
        outfile = open(output + '/' + output.split('/')[-1] + '_' + tf_val + '.sif',
                       'w')  # Generates sif output file for each TF
        outfile_tbl = open(output + '/' + output.split('/')[-1] + '_' + tf_val + '.tbl',
                           'w')  # Generates sif output file for
        # each TF
        sub_df[['TF', 'EDGE', 'TARGET']].to_csv(outfile, sep='\t', index=False, header=False)
        sub_df[['shared name', 'FOLDCHANGE', 'PVALUE']].to_csv(outfile_tbl, sep='\t', index=False)
        outfile.close()  # close the file resources
        outfile_tbl.close()

    total_exp = len(sif_rs_tabular.columns.tolist())  # count the total number of experiments

    sif_rs_tabular['target_count'] = (sif_rs_tabular.notnull() * 1).sum(axis=1)
    # df subset of common targets (i.e. targets repersented across all the exps)
    sub_common_targets = sif_rs_tabular[sif_rs_tabular['target_count'] == total_exp]
    # df subset of shared targets (i.e. target genes present in >1 exp)
    sub_shared_targets = sif_rs_tabular[sif_rs_tabular['target_count'] > 1]

    # I get warnings from pandas for using inplace with drop
    # sub_common_targets.drop('target_count', 1, inplace=True) # drop the target_count column
    # sub_shared_targets.drop('target_count', 1, inplace=True) # drop the target_count column
    # Alternative
    sub_common_targets = sub_common_targets.ix[:, 0:-1]  # removing the last columns target_count column
    sub_shared_targets = sub_shared_targets.ix[:, 0:-1]  # removing the last columns target_count column

    # Generates sif output file for all TF
    outfile_common = open(output + '/' + output.split('/')[-1] + '_commonTargets.sif', 'w')
    outfile_common_tbl = open(output + '/' + output.split('/')[-1] + '_commonTargets.tbl', 'w')
    stacked_common_targets = pd.DataFrame(sub_common_targets.stack().reset_index())
    stacked_common_targets.columns = ['TARGET', 'TF', 'EDGE']
    stacked_common_targets['TF'] = stacked_common_targets['TF'].apply(
        lambda x: x.split('_')[0])  # extract TFname from expID
    stacked_common_targets.drop_duplicates(
        inplace=True)  # dropping duplicates to solve the problem of multiple representation of
    # same experiment (with multiple analysis ID) and chipseq multiple time-points
    reordered_common_targets = stacked_common_targets.assign(
        **stacked_common_targets.EDGE.str.extract(regex, expand=True).to_dict('list'))
    reordered_common_targets[['TF', 'EDGE', 'TARGET']].to_csv(outfile_common, sep='\t', index=False, header=False)

    if not reordered_common_targets.empty:
        reordered_common_targets['shared name'] = reordered_common_targets['TARGET'] + ' (' + \
                                                  reordered_common_targets['EDGE'] + ') ' + \
                                                  reordered_common_targets['TF']
        reordered_common_targets[['shared name', 'FOLDCHANGE', 'PVALUE']].to_csv(outfile_common_tbl, sep='\t',
                                                                                 index=False)

    # Generates sif output file for shared targets
    outfile_shared = open(output + '/' + output.split('/')[-1] + '_sharedTargets.sif', 'w')
    outfile_shared_tbl = open(output + '/' + output.split('/')[-1] + '_sharedTargets.tbl', 'w')
    stacked_shared_targets = pd.DataFrame(sub_shared_targets.stack().reset_index())
    stacked_shared_targets.columns = ['TARGET', 'TF', 'EDGE']
    stacked_shared_targets['TF'] = stacked_shared_targets['TF'].apply(
        lambda x: x.split('_')[0])  # extract TFname from expID
    stacked_shared_targets.drop_duplicates(
        inplace=True)  # dropping duplicates to solve the problem of multiple representation of
    # same experiment (with multiple analysis ID) and chipseq multiple time-points
    reordered_shared_targets = stacked_shared_targets.assign(
        **stacked_shared_targets.EDGE.str.extract(regex, expand=True).to_dict('list'))
    reordered_shared_targets[['TF', 'EDGE', 'TARGET']].to_csv(outfile_shared, sep='\t', index=False)
    reordered_shared_targets['shared name'] = reordered_shared_targets['TARGET'] + ' (' + \
                                              reordered_shared_targets['EDGE'] + ') ' + \
                                              reordered_shared_targets['TF']
    reordered_shared_targets[['shared name', 'FOLDCHANGE', 'PVALUE']].to_csv(outfile_shared_tbl, sep='\t',
                                                                             index=False)

    outfile_all.close()  # close the file resources
    outfile_all_tbl.close()
    outfile_common.close()  # close the file resources
    outfile_common_tbl.close()
    outfile_shared.close()  # close the file resources
    outfile_shared_tbl.close()

    zip_name = shutil.make_archive(output, 'zip', output)  # create a zip file for output directory
    shutil.rmtree(output)  # delete the output directory after creating zip file

    return zip_name


#################################
# Generate genelists
def create_genelists(reordered_tmp_df, output):
    # read the pickled dataframe
    reordered_subdf= reordered_tmp_df[['TF','TARGET']]

    reordered_subdf_concat= (reordered_subdf.groupby('TF', as_index=False).apply(lambda x: pd.concat(
        [pd.Series('>' + x.iloc[0]['TF']),pd.Series(x.loc[:, 'TARGET'])])).reset_index(drop=True))
    #print('reordered_subdf_concat= ',reordered_subdf_concat)
    outfile_genelist = open(output + '/' + output.split('/')[-1] + '_genelist_allTFs.fa','w') # Generates genelists for all TF
    reordered_subdf_concat.to_csv(outfile_genelist, sep='\t', index=False, header=False)
    outfile_genelist.close()
