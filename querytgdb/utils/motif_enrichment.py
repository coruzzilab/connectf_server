import gzip
import pathlib
import pickle
import sys
from collections import OrderedDict
from functools import partial, reduce
from io import BytesIO
from itertools import chain, count, tee
from multiprocessing.pool import ThreadPool
from operator import or_
from typing import Dict, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as hierarchy
import seaborn as sns
from django.conf import settings
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import fdrcorrection

from querytgdb.models import Analysis
from ..utils import NAME_REGEX, column_string, svg_font_adder
from ..utils.parser import ANNOTATIONS


class MotifData:
    def __init__(self, data_dir: Union[str, pathlib.Path]):
        self.pool = ThreadPool()
        self._annotated_async = self.pool.apply_async(pd.read_pickle,
                                                      args=(data_dir / 'static' / 'annotated.pickle.gz',))
        self._annotated = None
        self._annotated_promo = None
        self._ann_promo_dedup = None
        self._promo_cluster_size = None

        self._annotated_body = None
        self._ann_body_dedup = None
        self._body_cluster_size = None

    def close(self):
        self.pool.terminate()

    @property
    def annotated(self):
        if self._annotated is None:
            annotated = self._annotated_async.get()
            self.close()
            annotated = annotated[annotated['p-value'] < 0.0001]

            self._annotated = annotated

        return self._annotated

    @property
    def annotated_promo(self):
        if self._annotated_promo is None:
            self._annotated_promo = self.annotated[
                (self.annotated['stop'] - self.annotated['start'] + self.annotated['dist']) < 0]

        return self._annotated_promo

    @property
    def ann_promo_dedup(self):
        if self._ann_promo_dedup is None:
            self._ann_promo_dedup = self.annotated_promo.drop_duplicates('match_id')

        return self._ann_promo_dedup

    @property
    def promo_cluster_size(self):
        if self._promo_cluster_size is None:
            self._promo_cluster_size = self.ann_promo_dedup.groupby('#pattern name').size()

        return self._promo_cluster_size

    @property
    def annotated_body(self):
        if self._annotated_body is None:
            self._annotated_body = self.annotated[self.annotated['dist'] > 0]

        return self._annotated_body

    @property
    def ann_body_dedup(self):
        if self._ann_body_dedup is None:
            self._ann_body_dedup = self.annotated_body.drop_duplicates('match_id')

        return self._ann_body_dedup

    @property
    def body_cluster_size(self):
        if self._body_cluster_size is None:
            self._body_cluster_size = self.ann_body_dedup.groupby('#pattern name').size()

        return self._body_cluster_size

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()


APPS_DIR = pathlib.Path(settings.BASE_DIR) / 'tgdbbackend'

MOTIF = MotifData(APPS_DIR)

CLUSTER_INFO = pd.read_csv(
    APPS_DIR / 'static' / 'cluster_info.csv.gz',
    index_col=0
).fillna('').to_dict('index')

COLORS = sns.color_palette("husl", 2)


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
                                                  annotated=MOTIF.annotated_promo,
                                                  annotated_dedup=MOTIF.ann_promo_dedup,
                                                  ann_cluster_size=MOTIF.promo_cluster_size), res.values()))
    if body:
        body_enrich, body_reject = zip(*map(partial(get_list_enrichment,
                                                    alpha=alpha,
                                                    annotated=MOTIF.annotated_body,
                                                    annotated_dedup=MOTIF.ann_body_dedup,
                                                    ann_cluster_size=MOTIF.body_cluster_size), res.values()))

        df = pd.concat(chain.from_iterable(zip(promo_enrich, body_enrich)), axis=1)
        _promo, _body = tee(res.keys(), 2)
        columns = list(chain.from_iterable(zip(
            map(lambda c: c + ('promo',), _promo),
            map(lambda c: c + ('body',), _body))))

        if show_reject:
            rejects = reduce(or_, chain(promo_reject, body_reject))
        else:
            rejects = pd.concat(chain.from_iterable(zip(promo_reject, body_reject)), axis=1)
    else:
        df = pd.concat(promo_enrich, axis=1, sort=True)
        columns = [c + ('promo',) for c in res.keys()]

        if show_reject:
            rejects = reduce(or_, promo_reject)
        else:
            rejects = pd.concat(promo_reject, axis=1, sort=True)

    df.columns = columns

    if show_reject:
        rejects = rejects.reindex(df.index)
        df = df[rejects]
    else:
        rejects.columns = columns
        df = df.where(rejects).dropna(how='all')

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

    tfs, exp_ids, analysis_ids = zip(*res.keys())
    analyses = Analysis.objects.filter(
        name__in=analysis_ids,
        experiment__name__in=exp_ids).prefetch_related('experiment')

    meta_dicts: OrderedDict[Tuple[str, ...], OrderedDict] = OrderedDict()

    for tf, exp_id, analysis_id in res.keys():
        m = NAME_REGEX.match(tf)

        if m:
            name, criterion = m.groups('')
        else:
            name = tf
            criterion = ''

        data = OrderedDict([('name', name), ('filter', criterion)])
        try:
            analysis = analyses.get(
                name=analysis_id,
                experiment__name=exp_id)
            data.update(analysis.analysisdata_set.values_list('key', 'value'))
            data.update(analysis.experiment.experimentdata_set.values_list('key', 'value'))
        except Analysis.DoesNotExist:
            pass
        meta_dicts[(tf, exp_id, analysis_id)] = data

    try:
        if target_genes_path:
            with gzip.open(target_genes_path, 'rb') as f:
                _, target_lists = pickle.load(f)
                meta_dicts.update(
                    [((name, *m_name), OrderedDict({'list_name': name, **m})) for m_name, m in meta_dicts.items() for
                     name in target_lists.keys()])
                # Use list comprehension instead of generator expression to avoid mutating dict as we update it
                res.update([((t_name, *r_name), t_list & r_list)
                            for r_name, r_list in res.items()
                            for t_name, t_list in target_lists.items()])
    except FileNotFoundError:
        pass

    df = motif_enrichment(res, alpha=alpha, show_reject=False, body=body)

    if df.empty:
        raise NoEnrichedMotif

    # rows, cols = df.shape
    # if rows > 1:
    #     z = hierarchy.linkage(df.values, method='average', optimal_ordering=True)
    #     df = df.iloc[hierarchy.leaves_list(z), :]
    #
    # if cols > 1:
    #     z = hierarchy.linkage(df.values.T, method='average', optimal_ordering=True)
    #     df = df.iloc[:, hierarchy.leaves_list(z)]

    df = df.where(pd.notnull(df), None)

    return {
        'columns': [meta_dicts[c] for c in res.keys()],
        'result': list(merge_cluster_info(df))
    }


def get_motif_enrichment_heatmap(cache_path, target_genes_path=None, alpha=0.05, lower_bound=None, upper_bound=None,
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
    df = df.clip_upper(sys.maxsize)

    df = df.rename(index={idx: "{} ({})".format(idx, CLUSTER_INFO[idx]['Family']) for idx in df.index})

    names, exp_ids, analysis_ids = zip(*res.keys())

    analyses = Analysis.objects.filter(
        name__in=analysis_ids,
        experiment__name__in=exp_ids).prefetch_related('experiment')

    columns = []

    for (name, exp_id, analysis_id), col_name in zip(res.keys(), map(column_string, count(1))):
        tf = analyses.get(name=analysis_id, experiment__name=exp_id).experiment.tf
        columns.append('{1} — {0}{2}'.format(tf.gene_id, col_name, f' ({tf.name})' if tf.name else ''))

    if not body:
        df.columns = columns
    else:
        _promo, _body = tee(columns, 2)
        df.columns = chain.from_iterable(zip(
            map(lambda x: x + ' promo', _promo),
            map(lambda x: x + ' body', _body)
        ))

    # df = df.rename(columns='_'.join)

    df = df.T

    rows, cols = df.shape

    opts = {}
    if rows > 1:
        opts['row_linkage'] = hierarchy.linkage(df.values, method='average', optimal_ordering=True)
        if body:
            opts['row_colors'] = [COLORS[s.endswith('promo')] for s in df.index]

    if cols > 1:
        opts['col_linkage'] = hierarchy.linkage(df.values.T, method='average', optimal_ordering=True)

    plt.figure()
    heatmap_graph = sns.clustermap(df,
                                   cmap="YlGnBu",
                                   cbar_kws={'label': 'Enrichment (-log10 p)'},
                                   xticklabels=1,
                                   yticklabels=1,
                                   row_cluster=rows > 1,
                                   col_cluster=cols > 1,
                                   vmin=lower_bound,
                                   vmax=upper_bound,
                                   **opts)
    plt.setp(heatmap_graph.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    plt.setp(heatmap_graph.ax_heatmap.xaxis.get_majorticklabels(), rotation=270)

    buff = BytesIO()
    heatmap_graph.savefig(buff)
    buff.seek(0)
    svg_font_adder(buff)
    buff.seek(0)

    return buff


def get_motif_enrichment_heatmap_table(cache_path, target_genes_path=None):
    df = pd.read_pickle(cache_path)
    df = df.loc[:, (slice(None), slice(None), slice(None), 'EDGE')]
    df.columns = df.columns.droplevel(3)
    res = df.columns.tolist()

    try:
        if target_genes_path:
            with gzip.open(target_genes_path, 'rb') as f:
                _, target_lists = pickle.load(f)
                res.extend([(t_name, *r_name)
                            for r_name in res
                            for t_name in target_lists.keys()])
    except FileNotFoundError:
        pass

    names, exp_ids, analysis_ids = zip(*res)

    analyses = Analysis.objects.filter(
        name__in=analysis_ids,
        experiment__name__in=exp_ids).prefetch_related('experiment')

    for (name, exp_id, analysis_id), col_str in zip(res, map(column_string, count(1))):
        m = NAME_REGEX.match(name)

        if m:
            name, criterion = m.groups('')
        else:
            criterion = ''

        info = OrderedDict([('name', name), ('filter', criterion)])

        try:
            analysis = analyses.get(name=analysis_id, experiment__name=exp_id)
            gene_id = analysis.experiment.tf.gene_id

            try:
                gene_name = ANNOTATIONS.at[gene_id, 'Name']
            except KeyError:
                gene_name = ''

            info.update(analysis.experiment.experimentdata_set.values_list('key', 'value'))
            info.update(analysis.analysisdata_set.values_list('key', 'value'))

            yield (info, col_str, name, criterion, gene_name, analysis_id)
        except Analysis.DoesNotExist:
            pass
