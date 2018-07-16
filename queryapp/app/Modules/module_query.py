'''
Call this module to query the TargetDB database. Returns a pandas dataframe.

    My main program query_TargetDb calls the function queryTFDB to generate df
    for each TF
'''

##############
# Modules
import pandas as pd
from sqlalchemy.orm import aliased

from .create_mysqlDB import Edges, Interactions, Nodes


################################################
# Query the database
def queryTFDB(sess, q_TFname, qedgelist, qmetalist):
    # create alias for node1
    n1 = aliased(Nodes)
    n2 = aliased(Nodes)
    count = 0
    list_metaid = list()

    # print 'q_TFname= ',q_TFname
    # print 'qedgelist= ',qedgelist
    # print 'qmetalist= ',qmetalist

    if qedgelist:
        if qmetalist:  # if both edges and meta ids are given in the user query
            rs = sess.query(n1.node_name, Edges.edge_name, n2.node_name,
                            Interactions.meta_id,
                            Interactions.analysis_id).filter(
                n1.node_name == q_TFname). \
                filter(n1.node_id == Interactions.node_1_id).filter(
                Interactions.node_2_id == n2.node_id). \
                filter(Interactions.edge_id == Edges.edge_id).filter(
                Edges.edge_name.in_(qedgelist)).filter(
                Interactions.meta_id.in_(qmetalist)).all()
        if not qmetalist:  # if only edges are given in the user query
            rs = sess.query(n1.node_name, Edges.edge_name, n2.node_name,
                            Interactions.meta_id,
                            Interactions.analysis_id).filter(
                n1.node_name == q_TFname). \
                filter(n1.node_id == Interactions.node_1_id).filter(
                Interactions.node_2_id == n2.node_id). \
                filter(Interactions.edge_id == Edges.edge_id).filter(
                Edges.edge_name.in_(qedgelist)).all()
    if qmetalist and not qedgelist:  # if only metaid is given in the user query
        rs = sess.query(n1.node_name, Edges.edge_name, n2.node_name,
                        Interactions.meta_id, Interactions.analysis_id).filter(
            n1.node_name == q_TFname). \
            filter(n1.node_id == Interactions.node_1_id).filter(
            Interactions.node_2_id == n2.node_id). \
            filter(Interactions.edge_id == Edges.edge_id).filter(
            Interactions.meta_id.in_(qmetalist)).all()
    if not qmetalist and not qedgelist:  # if only TFs are given- no edges
        # and no metadata asked in the query
        rs = sess.query(n1.node_name, Edges.edge_name, n2.node_name,
                        Interactions.meta_id, Interactions.analysis_id).filter(
            n1.node_name == q_TFname). \
            filter(n1.node_id == Interactions.node_1_id).filter(
            Interactions.node_2_id == n2.node_id). \
            filter(Interactions.edge_id == Edges.edge_id).all()

    rs_pd = pd.DataFrame(rs,
                         columns=['TF', 'EDGE', 'TARGET', 'META', 'ANALYSIS'])
    #	if rs_pd.empty:
    #		print 'No data matched your query!'
    #		sys.exit(1)

    rs_pd['ANALYSIS'] = rs_pd['ANALYSIS'].str.replace('.',
                                                      '_')  # pvalues '.' are
    #  replaces because pandas does not allow to use these chars with
    # pandas.query
    rs_pd['META_ID'] = rs_pd[['META', 'ANALYSIS']].apply(lambda x: '_'.join(x),
                                                         axis=1)  #
    # separating metaid and analysis id with '_' because pandas does not
    # allow other chars
    rs_pd.drop('META', axis=1, inplace=True)
    rs_pd.drop('ANALYSIS', axis=1, inplace=True)

    if not rs_pd.empty:
        # code to change CHIPSEQ column names for time points: Example:
        # METAID_X will be METAID_X:0, METAID_X:1, METAID_X:5
        pattern = rs_pd.META_ID.str.contains('CHIPSEQ')
        # attach time point with each metaid
        rs_pd.loc[pattern, 'META_ID'] = rs_pd.loc[pattern, 'META_ID'] + '_' + \
                                        rs_pd.loc[pattern, 'EDGE'].str.split(
                                            ':').str.get(2)
    return rs_pd
