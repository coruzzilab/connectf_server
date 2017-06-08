'''
Call this module when ALLRNASEQ, Allinduced, Allrepressed and Allchipseq
'''

def create_edges_allquery():
    all_edge_flag = all_edge_flag + 1
    orig_edges = tmp_edge.split()

    for pos, val_q_ed in enumerate(edges):
        allquery = None
        if 'ALLRNASEQ' in val_q_ed.upper():  # Get ALLRNASEQ
            # edges. This will not work if AND [ALLRNASEQ] given
            # for edges. Go for ALLINDUCED or ALLREPRESSED in
            # that case.
            discard_pos.append(str(pos - 1))
            allquery = '%RNASEQ%'
            edges = create_query(sess, edges, allquery, q_tf,
                                     pos)
        if 'ALLCHIPSEQ' in val_q_ed.upper():  # Get ALLCHIPSEQ
            # edges.
            discard_pos.append(str(pos - 1))
            allquery = '%CHIPSEQ%'
            edges = create_query(sess, edges, allquery, q_tf,
                                     pos)
        if 'ALLINDUCED' in val_q_ed.upper():  # Get ALLINDUCED
            # edges. This will fetch edges for pCHX as well as
            # mCHX (if present for a TF).
            discard_pos.append(str(pos - 1))
            allquery = '%INDUCED%'
            edges = create_query(sess, edges, allquery, q_tf,
                                     pos)
        if 'ALLREPRESSED' in val_q_ed.upper():  # Get
            # ALLREPRESSED edges i.e. Edges for pCHX as well as
            # mCHX.
            discard_pos.append(str(pos - 1))
            allquery = '%REPRESSED%'
            edges = create_query(sess, edges, allquery, q_tf,
                                     pos)

    for p in sorted(discard_pos,reverse=True):  # this will fail in case when
        #  both chipseq and rnaseq are provided- TEST
        del edges[int(
            p)]  # delete all the positions with 'and' 'or' in
        # reverse order to avoid change in the pos after deleting

    return edges


###################################################################################
# Function to create queries when allrnaseq, allchipseq, allinduced or
# allrepressed is given in the user query
def create_query(sess, edges, allquery, q_tf, pos):
    # This code query the database and fetches edges for a given TF (q_TF)
    # and replaces allRNAseq (or allchipseq) with or/and between set of edges
    #  for a TF
    no1 = aliased(Nodes)
    no2 = aliased(Nodes)

    allquery_edges = list()
    tmpresult_edges = set(
        sess.query(Edges.edge_name).filter(no1.node_name == q_tf). \
            filter(no1.node_id == Interactions.node_1_id).filter(
            Interactions.node_2_id == no2.node_id). \
            filter(Interactions.edge_id == Edges.edge_id).filter(
            Edges.edge_name.like(allquery)).all())

    for x_r in tmpresult_edges:
        allquery_edges.append(x_r[0])
    all_r_edges = (' ' + edges[pos - 1].strip().replace('[', '') + ' ').join(
        allquery_edges)
    edges[pos] = '[' + all_r_edges + ']'
    return edges  # return replaced edge