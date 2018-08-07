import re
from itertools import groupby
from operator import itemgetter
from typing import List, Sized

import numpy as np
import pandas as pd

from ..utils.parser import expand_ref_ids

__all__ = ('create_export_zip',)


def column_string(n):
    s = ""
    while n > 0:
        n, remainder = divmod(n - 1, 26)
        s = chr(65 + remainder) + s
    return s


def create_export_zip(cache_dir, out_dir):
    df = pd.read_pickle(cache_dir + '/tabular_output.pickle.gz')
    create_sifs(df, out_dir)
    create_all_tf_genelists(df, out_dir)
    write_excel(cache_dir, out_dir)


def write_excel(cache_dir, out_dir):
    with pd.ExcelWriter(out_dir + '/tabular_output.xlsx', engine='xlsxwriter') as writer:
        df = pd.read_pickle(cache_dir + '/formatted_tabular_output.pickle.gz')
        metadata = pd.read_pickle(cache_dir + '/metadata.pickle.gz')
        metadata.index.name = None

        write_data(df, writer)
        write_metadata(metadata, writer)


def write_metadata(df: pd.DataFrame, writer: pd.ExcelWriter):
    df = expand_ref_ids(df)
    df.to_excel(writer, 'metadata', header=True)

    word_lens = df.astype(str).applymap(len)
    max_len = np.nanmax(word_lens[(word_lens - word_lens.values.mean()).abs() < 2 * word_lens.values.std()].values)

    workbook = writer.book
    worksheet = writer.sheets['metadata']

    header_fmt = workbook.add_format(
        {'font_name': 'Calibri', 'font_size': 15, 'bold': True, 'align': 'center', 'border': 1})

    worksheet.set_row(1, None, header_fmt)
    worksheet.set_column(0, df.shape[1], max_len)


def write_data(data: List[Sized], writer: pd.ExcelWriter):
    workbook = writer.book
    worksheet = workbook.add_worksheet('table')

    row_num, col_num = len(data), len(data[0])

    bold_font = workbook.add_format({'bold': True, 'font_size': 13, 'border': 1, 'align': 'center'})

    for i, row in enumerate(data[6:]):
        worksheet.write_row(i + 6, 0, row)

    worksheet.set_column(0, col_num - 1, 15)
    worksheet.set_column(6, 7, None, bold_font)
    worksheet.set_column(6, 6, 10)

    last_col = column_string(col_num)

    format_repressed = workbook.add_format({'bg_color': '#FA8072', 'font_color': '#000000'})
    format_induced = workbook.add_format({'bg_color': '#98FB98', 'font_color': '#000000'})
    format_present = workbook.add_format({'bg_color': '#F7DC6F', 'font_color': '#000000'})
    merge_format = workbook.add_format({
        'bold': 1,
        'font_size': 14,
        'border': 1,
        'align': 'center',
        'valign': 'vcenter',
        'fg_color': 'white'})

    worksheet.conditional_format('I7:{}{}'.format(last_col, row_num),
                                 {'type': 'text', 'criteria': 'containing', 'value': 'INDUCED',
                                  'format': format_induced})
    worksheet.conditional_format('I7:{}{}'.format(last_col, row_num),
                                 {'type': 'text', 'criteria': 'containing', 'value': 'REPRESSED',
                                  'format': format_repressed})
    worksheet.conditional_format('I7:{}{}'.format(last_col, row_num),
                                 {'type': 'text', 'criteria': 'containing', 'value': 'Present',
                                  'format': format_present})

    columns = list(zip(*data[:6]))
    for i in range(6):
        index = 0
        for label, group in groupby(columns, key=itemgetter(*range(0, i + 1))):
            size = sum(1 for _ in group)

            if isinstance(label, str):
                s = label
            else:
                s = label[-1]

            if size > 1:
                worksheet.merge_range(i, index, i, index + size - 1, s, merge_format)
            else:
                worksheet.write(i, index, s, merge_format)

            index += size


group_opt = {
    'axis': 1,
    'level': 0,
    'by': lambda x: re.sub(r'_\w+$', '', x)
}


#################################
# Generate sif output
def create_sifs(result, output):
    res_group = result.groupby(**group_opt)

    with open(output + '/all_tfs.sif', 'w') as all_sif, open(output + '/all_tfs.tbl', 'w') as all_tbl:
        all_tbl.write('shared name\tFOLDCHANGE\tPVALUE\n')

        for name, group in res_group:
            g = group.stack(level=[0, 1, 2]).reset_index(level=[1, 2, 3], drop=True).reset_index().drop_duplicates()

            with open(output + '/{}.sif'.format(name), 'w') as f:
                temp_sif = g.groupby('EDGE').apply(
                    lambda x: '{}\t{}\t{}'.format(name, x.name, '\t'.join(x['TARGET']))
                )
                temp_sif.to_csv(f, index=False)
                temp_sif.to_csv(all_sif, index=False)
            with open(output + '/{}.tbl'.format(name), 'w') as f:
                f.write('shared name\tFOLDCHANGE\tPVALUE\n')
                temp_tbl = g.apply(lambda x: '{0} ({1[EDGE]}) {1[TARGET]}\t{1[Log2FC]}\t{1[Pvalue]}'.format(name, x),
                                   axis=1)
                temp_tbl.to_csv(f, index=False)
                temp_tbl.to_csv(all_tbl, index=False)

    crit_tfs = res_group.apply(lambda x: x.loc[:, (slice(None), slice(None), slice(None), 'EDGE')].notna().any(axis=1))
    common_tfs = crit_tfs.all(axis=1)
    shared_tfs = crit_tfs.sum(axis=1) > 1

    def create_filtered_sifs(index, sif_file, tbl_file):
        with open(sif_file, 'w') as sif, open(tbl_file, 'w') as tbl:
            tbl.write('shared name\tFOLDCHANGE\tPVALUE\n')
            for name, group in result.loc[index, :].groupby(**group_opt):
                g = group.stack(level=[0, 1, 2]).reset_index(level=[1, 2, 3], drop=True).reset_index().drop_duplicates()
                # with pd.option_context('display.width', 200, 'display.max_columns', 20, 'display.max_rows', 20000):

                temp_sif = g.groupby('EDGE').apply(
                    lambda x: '{}\t{}\t{}'.format(name, x.name, '\t'.join(x['TARGET']))
                )
                temp_sif.to_csv(sif, index=False)
                temp_tbl = g.apply(lambda x: '{0} ({1[EDGE]}) {1[TARGET]}\t{1[Log2FC]}\t{1[Pvalue]}'.format(name, x),
                                   axis=1)
                temp_tbl.to_csv(tbl, index=False)

    create_filtered_sifs(common_tfs, output + '/common_tfs.sif', output + '/common_tfs.tbl')
    create_filtered_sifs(shared_tfs, output + '/shared_tfs.sif', output + '/shared_tfs.tbl')


#################################
# Generate genelists
def create_all_tf_genelists(result: pd.DataFrame, output):
    group_res = result.groupby(**group_opt)

    with open(output + '/genelists_all_tf.fa', 'w') as f:
        for name, group in group_res:
            f.write('>{}\n{}\n'.format(name, '\n'.join(group.index)))
