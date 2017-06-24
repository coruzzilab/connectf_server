'''
Call this module to query the TargetDB database. Returns a pandas dataframe.

    My main program query_TargetDb calls the function queryTFDB to generate df
    for each TF
'''

##############
# Modules
from __future__ import absolute_import

import pandas as pd

from ....models import Interactions, ReferenceId


# from create_mysqlDB_v2 import TargetDBTF, Edges, Annotation, Interactions, ReferenceId, Metadata, Analysis

################################################
# Query the database
# @profile
def queryTFDB(q_TFname, rs_meta_list):
    rs = list(Interactions.objects.select_related().filter(db_tf_id__db_tf_agi__exact=q_TFname).values_list(
        'db_tf_id__db_tf_agi', 'edge_id__edge_name', 'target_id__agi_id', 'ref_id__ref_id'))

    rs_pd_tmp = pd.DataFrame(rs, columns=['TF', 'EDGE', 'TARGET', 'REFID'])
    if rs_meta_list:
        rs_pd = rs_pd_tmp.loc[rs_pd_tmp.REFID.isin(rs_meta_list)]
    else:
        rs_pd = rs_pd_tmp

    list_ref_id = rs_pd.REFID.unique()
    # print('list_ref_id= ',list_ref_id)
    meta_ref = ReferenceId.objects.select_related().filter(ref_id__in=list_ref_id). \
        values_list('ref_id', 'meta_id__meta_fullid', 'analysis_id__analysis_fullid')
    # print('meta_ref= ',meta_ref)
    # print('meta_ref= ', meta_ref)

    meta_ref_dict = dict()
    for val_m in meta_ref:
        meta_ref_dict[val_m[0]] = '_'.join([val_m[1], val_m[2], str(val_m[0])])

    # print('meta_ref_dict= ',meta_ref_dict)

    # Pandas query func throws an error if columns names are numbers so I had to include meta_id in RefID
    # column name '1', '1_2', '1_a' etc. will not work
    if not rs_pd.empty:
        rs_pd.REFID.replace(to_replace=meta_ref_dict, inplace=True)
        # pvalues '.' are replaces because pandas does not allow to use these chars with pandas.query
        rs_pd['REFID'] = rs_pd['REFID'].str.replace('.', '_')

    # print('rs_pd= ',rs_pd)
    return rs_pd
