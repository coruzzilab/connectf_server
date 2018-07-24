import gzip
import pathlib
import pickle
from collections import OrderedDict
from functools import partial, reduce
from io import BytesIO
from itertools import chain
from operator import or_
from typing import Dict, Tuple

import matplotlib
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as hierarchy
from django.conf import settings
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import fdrcorrection

from querytgdb.models import ReferenceId

matplotlib.use('SVG')
import matplotlib.pyplot as plt
import seaborn as sns

APPS_DIR = pathlib.Path(settings.BASE_DIR) / 'tgdbbackend'

ANNOTATED = pd.read_pickle(APPS_DIR / 'static' / 'annotated.pickle.gz')
ANNOTATED = ANNOTATED[ANNOTATED['p-value'] < 0.0001]

ANNOTATED_PROMO = ANNOTATED[(ANNOTATED['stop'] - ANNOTATED['start'] + ANNOTATED['dist']) < 0]
ANN_PROMO_DEDUP = ANNOTATED_PROMO.drop_duplicates('match_id')
PROMO_CLUSTER_SIZE = ANN_PROMO_DEDUP.groupby('#pattern name').size()

ANNOTATED_BODY = ANNOTATED[ANNOTATED['dist'] > 0]
ANN_BODY_DEDUP = ANNOTATED_BODY.drop_duplicates('match_id')
BODY_CLUSTER_SIZE = ANN_BODY_DEDUP.groupby('#pattern name').size()

CLUSTER_INFO = pd.read_pickle(
    APPS_DIR / 'static' / 'cluster_info.pickle.gz',
    compression='gzip'
).to_dict('index')


class NoEnrichedMotif(ValueError):
    pass


def motif_enrichment(res: Dict[Tuple[str], pd.Series], alpha: float = 0.05, show_reject: bool = True,
                     body: bool = False) -> pd.DataFrame:
    def get_list_enrichment(gene_list, annotated, annotated_dedup, ann_cluster_size,
                            alpha: float = 0.05) -> Tuple[pd.Series, pd.Series]:
        list_cluster_dedup = annotated[annotated.index.isin(gene_list)].drop_duplicates('match_id')
        list_cluster_size = list_cluster_dedup.groupby('#pattern name').size()

        def cluster_fisher(row):
            return fisher_exact(
                [[row[0], row[1] - row[0]],
                 [list_cluster_dedup.shape[0] - row[0],
                  annotated_dedup.shape[0] - list_cluster_dedup.shape[0] - row[1] + row[0]]],
                alternative='greater')[1]

        p_values = pd.concat([list_cluster_size, ann_cluster_size],
                             axis=1, sort=False).fillna(0).apply(cluster_fisher, axis=1).sort_values()
        reject, adj_p = fdrcorrection(p_values, alpha=alpha, is_sorted=True)

        str_index = p_values.index.astype(str)

        return pd.Series(adj_p, index=str_index), pd.Series(reject, index=str_index)

    promo_enrich, promo_reject = zip(*map(partial(get_list_enrichment,
                                                  alpha=alpha,
                                                  annotated=ANNOTATED_PROMO,
                                                  annotated_dedup=ANN_PROMO_DEDUP,
                                                  ann_cluster_size=PROMO_CLUSTER_SIZE), res.values()))
    if body:
        body_enrich, body_reject = zip(*map(partial(get_list_enrichment,
                                                    alpha=alpha,
                                                    annotated=ANNOTATED_BODY,
                                                    annotated_dedup=ANN_BODY_DEDUP,
                                                    ann_cluster_size=BODY_CLUSTER_SIZE), res.values()))

        df = pd.concat(chain.from_iterable(zip(promo_enrich, body_enrich)), axis=1)
        columns = list(chain.from_iterable(zip(
            map(lambda c: '_'.join(c) + '_promo', res.keys()),
            map(lambda c: '_'.join(c) + '_body', res.keys()))))

        if show_reject:
            rejects = reduce(or_, chain(promo_reject, body_reject))
        else:
            rejects = pd.concat(chain.from_iterable(zip(promo_reject, body_reject)), axis=1)
    else:
        df = pd.concat(promo_enrich, axis=1, sort=True)
        columns = ['_'.join(c) + '_promo' for c in res.keys()]

        if show_reject:
            rejects = reduce(or_, promo_reject)
        else:
            rejects = pd.concat(promo_reject, axis=1, sort=True)

    df.columns = columns

    if show_reject:
        df = df[rejects]
    else:
        rejects.columns = columns
        df[~rejects] = np.nan
        df.dropna(how='all', inplace=True)

    return df


def merge_cluster_info(df):
    for idx, *row in df.itertuples(name=None):
        info = {'name': idx}
        try:
            info.update(CLUSTER_INFO[idx])
        except KeyError:
            pass
        yield [info] + row


def get_motif_enrichment_json(cache_path, target_genes_path=None, alpha=0.05, body=False) -> Dict:
    df = pd.read_pickle(cache_path)
    # remove unneeded info from dataframe
    df = df.loc[:, (slice(None), slice(None), slice(None), 'EDGE')]
    df.columns = df.columns.droplevel(3)
    res = OrderedDict((name, set(col.index[col.notnull()])) for name, col in df.iteritems())

    tfs, meta_ids, analysis_ids = zip(*res.keys())
    references = ReferenceId.objects.filter(analysis_id__analysis_fullid__in=analysis_ids,
                                            meta_id__meta_fullid__in=meta_ids)

    meta_dicts = []

    for meta_id, analysis_id in zip(meta_ids, analysis_ids):
        data = OrderedDict()
        try:
            ref = references.get(
                analysis_id__analysis_fullid=analysis_id,
                meta_id__meta_fullid=meta_id)
            data.update(OrderedDict(
                ref.analysis_id.analysisiddata_set.values_list('analysis_type', 'analysis_value')))
            data.update(
                ref.meta_id.metaiddata_set.values_list('meta_type', 'meta_value'))
        except ReferenceId.DoesNotExist:
            pass
        meta_dicts.append(data)

    try:
        if target_genes_path:
            with gzip.open(target_genes_path, 'rb') as f:
                _, target_lists = pickle.load(f)
                meta_dicts.extend(
                    [{'list_name': name, **m} for m in meta_dicts for name in target_lists.keys()])
                # Use list comprehension instead of generator expression to avoid mutating dict as we update it
                res.update([((t_name, *r_name), t_list & r_list)
                            for r_name, r_list in res.items()
                            for t_name, t_list in target_lists.items()])
    except FileNotFoundError:
        pass

    df = motif_enrichment(res, alpha=alpha, show_reject=False, body=body)

    if df.empty:
        raise NoEnrichedMotif

    df = df.where(pd.notnull(df), None)

    return {
        'columns': OrderedDict(zip(('_'.join(r) for r in res.keys()), meta_dicts)),
        'result': list(merge_cluster_info(df))
    }


def get_motif_enrichment_heatmap(cache_path, target_genes_path=None, alpha=0.05, lower_bound=0, upper_bound=10,
                                 body=False) -> BytesIO:
    df = pd.read_pickle(cache_path)
    df = df.loc[:, (slice(None), slice(None), slice(None), 'EDGE')]
    df.columns = df.columns.droplevel(3)
    res = OrderedDict((name, set(col.index[col.notnull()])) for name, col in df.iteritems())

    try:
        if target_genes_path:
            with gzip.open(target_genes_path, 'rb') as f:
                _, target_lists = pickle.load(f)
                res.update([((t_name, *r_name), t_list & r_list)
                            for r_name, r_list in res.items()
                            for t_name, t_list in target_lists.items()])
    except FileNotFoundError:
        pass

    df = motif_enrichment(res, alpha=alpha, body=body)
    if df.empty:
        raise NoEnrichedMotif
    df = -np.log10(df)

    if lower_bound:
        df[df < lower_bound] = lower_bound
    # make max 10 for overly small p-values
    df[df > upper_bound] = upper_bound

    df = df.rename(index={idx: "{} ({})".format(idx, CLUSTER_INFO[idx]['Family']) for idx in df.index})

    rows, cols = df.shape

    opts = {}
    if rows > 1:
        opts['row_linkage'] = hierarchy.linkage(df.values, method='average', optimal_ordering=True)
    if cols > 1:
        opts['col_linkage'] = hierarchy.linkage(df.values.T, method='average', optimal_ordering=True)

    heatmap_graph = sns.clustermap(df,
                                   cmap="YlGnBu",
                                   xticklabels=1,
                                   row_cluster=rows > 1,
                                   col_cluster=cols > 1,
                                   **opts)
    plt.setp(heatmap_graph.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    plt.setp(heatmap_graph.ax_heatmap.xaxis.get_majorticklabels(), rotation=270)

    buff = BytesIO()
    heatmap_graph.savefig(buff)
    buff.seek(0)

    return buff
