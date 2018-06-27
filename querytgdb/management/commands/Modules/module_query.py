"""
Call this module to query the TargetDB database. Returns a pandas dataframe.

    My main program query_TargetDb calls the function queryTFDB to generate df
    for each TF
"""
from __future__ import absolute_import

import pandas as pd

from ....models import Interactions, ReferenceId


# from create_mysqlDB_v2 import TargetDBTF, Edges, Annotation, Interactions, ReferenceId, Metadata, Analysis

################################################
# Query the database
# @profile
def queryTFDB(q_TFname, rs_meta_list):
    rs_pd = pd.DataFrame(
        Interactions.objects.select_related().filter(db_tf_id__db_tf_agi__exact=q_TFname).values_list(
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
