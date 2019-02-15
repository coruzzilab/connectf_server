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
from typing import Dict, Generator, Iterable, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as hierarchy
import seaborn as sns
from django.conf import settings
from django.core.exceptions import ObjectDoesNotExist
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import fdrcorrection

from querytgdb.models import Analysis
from ..utils import clear_data, column_string, split_name, svg_font_adder
from ..utils.parser import annotations


class MotifData:
    def __init__(self, data_file: Union[str, pathlib.Path]):
        self.pool = ThreadPool()
        self._annotated_async = self.pool.apply_async(pd.read_pickle, args=(data_file,))
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


MOTIF = MotifData(settings.MOTIF_ANNOTATION)

CLUSTER_INFO = pd.read_csv(
    settings.MOTIF_CLUSTER,
    index_col=0
).fillna('').to_dict('index')

COLORS = sns.color_palette("husl", 2)


class NoEnrichedMotif(ValueError):
    pass


def cluster_fisher(row):
    return fisher_exact((row[:2], row[2:]), alternative='greater')[1]


def get_list_enrichment(gene_list, annotated, annotated_dedup, ann_cluster_size,
                        alpha: float = 0.05) -> Tuple[pd.Series, pd.Series]:
    list_cluster_dedup = annotated[annotated.index.isin(gene_list)].drop_duplicates('match_id')
    list_cluster_size = list_cluster_dedup.groupby('#pattern name').size()

    p_values = pd.concat([
        list_cluster_size,
        ann_cluster_size - list_cluster_size,
        list_cluster_dedup.shape[0] - list_cluster_size,
        annotated_dedup.shape[0] - list_cluster_dedup.shape[0] - ann_cluster_size + list_cluster_size
    ], axis=1, sort=False).fillna(0).apply(cluster_fisher, axis=1).sort_values()

    reject, adj_p = fdrcorrection(p_values, alpha=alpha, is_sorted=True)

    str_index = p_values.index.astype(str)

    return pd.Series(adj_p, index=str_index), pd.Series(reject, index=str_index)


def motif_enrichment(res: Dict[Tuple[str, Union[None, int]], Iterable], alpha: float = 0.05, show_reject: bool = True,
                     body: bool = False) -> pd.DataFrame:
    promo_enrich, promo_reject = zip(*map(partial(get_list_enrichment,
                                                  alpha=alpha,
                                                  annotated=MOTIF.annotated_promo,
                                                  annotated_dedup=MOTIF.ann_promo_dedup,
                                                  ann_cluster_size=MOTIF.promo_cluster_size),
                                          res.values()))
    if body:
        body_enrich, body_reject = zip(*map(partial(get_list_enrichment,
                                                    alpha=alpha,
                                                    annotated=MOTIF.annotated_body,
                                                    annotated_dedup=MOTIF.ann_body_dedup,
                                                    ann_cluster_size=MOTIF.body_cluster_size),
                                            res.values()))

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
    df = clear_data(df)
    res = OrderedDict((name, set(col.index[col.notnull()])) for name, col in df.iteritems())

    try:
        with gzip.open(target_genes_path, 'rb') as f:
            _, target_lists = pickle.load(f)
            res[('Target Gene List', None)] = set(chain.from_iterable(target_lists.values()))
    except (FileNotFoundError, TypeError):
        pass

    tfs, analysis_ids = zip(*res.keys())
    analyses = Analysis.objects.prefetch_related('tf').filter(pk__in=analysis_ids)

    columns = []

    for (tf, analysis_id), value in res.items():
        name, criterion = split_name(tf)

        data = OrderedDict([('name', name), ('filter', criterion)])
        try:
            analysis = analyses.get(pk=analysis_id)
            data['gene_name'] = analysis.tf.name
            data.update(analysis.meta_dict)
        except ObjectDoesNotExist:
            if (tf, analysis_id) == ('Target Gene List', None):
                data['genes'] = ", ".join(value)
            else:
                raise

        columns.append(data)

    df = motif_enrichment(res, alpha=alpha, show_reject=False, body=body)

    if df.empty:
        raise NoEnrichedMotif

    df = df.where(pd.notnull(df), None)

    return {
        'columns': columns,
        'result': list(merge_cluster_info(df))
    }


def get_motif_enrichment_heatmap(cache_path, target_genes_path=None, alpha=0.05, lower_bound=None, upper_bound=None,
                                 body=False) -> BytesIO:
    df = pd.read_pickle(cache_path)
    df = clear_data(df)
    res = OrderedDict((name, set(col.index[col.notnull()])) for name, col in df.iteritems())

    try:
        with gzip.open(target_genes_path, 'rb') as f:
            _, target_lists = pickle.load(f)
            res[('Target Gene List', None)] = set(chain.from_iterable(target_lists.values()))
    except (FileNotFoundError, TypeError):
        pass

    df = motif_enrichment(res, alpha=alpha, body=body)

    if df.empty:
        raise NoEnrichedMotif

    df = df.clip_lower(sys.float_info.min)
    df = -np.log10(df)

    df = df.rename(index={idx: "{} ({})".format(idx, CLUSTER_INFO[idx]['Family']) for idx in df.index})

    names, analysis_ids = zip(*res.keys())

    analyses = Analysis.objects.filter(pk__in=analysis_ids)

    columns = []

    for (name, analysis_id), col_name in zip(res.keys(), map(column_string, count(1))):
        try:
            tf = analyses.get(pk=analysis_id).tf
            columns.append('{1} — {0}{2}'.format(tf.gene_id, col_name, f' ({tf.name})' if tf.name else ''))
        except Analysis.DoesNotExist:
            columns.append('{} — {}'.format(col_name, name))

    if not body:
        df.columns = columns
    else:
        _promo, _body = tee(columns, 2)
        df.columns = chain.from_iterable(zip(
            map(lambda x: x + ' promo', _promo),
            map(lambda x: x + ' body', _body)
        ))

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
                                   method="ward",
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
    plt.close()
    buff.seek(0)
    svg_font_adder(buff)

    return buff


TableRow = Tuple[Dict, str, str, str, str, str]


def get_motif_enrichment_heatmap_table(cache_path, target_genes_path=None) -> Generator[TableRow, None, None]:
    df = pd.read_pickle(cache_path)
    df = clear_data(df)
    res = df.columns.tolist()

    names, analysis_ids = zip(*res)

    analyses = Analysis.objects.filter(pk__in=analysis_ids)

    col_strings = map(column_string, count(1))  # this is an infinite generator lazy evaluation only

    for (name, analysis_id), col_str in zip(res, col_strings):
        name, criterion = split_name(name)

        info = OrderedDict([('name', name), ('filter', criterion)])

        try:
            analysis = analyses.get(pk=analysis_id)
            gene_id = analysis.tf.gene_id

            try:
                gene_name = annotations().at[gene_id, 'Name']
            except KeyError:
                gene_name = ''

            info.update(analysis.analysisdata_set.values_list('key__name', 'value'))
        except Analysis.DoesNotExist:
            gene_name = ''

        yield (info, col_str, name, criterion, gene_name, analysis_id)

    try:
        with gzip.open(target_genes_path, 'rb') as f:
            _, target_lists = pickle.load(f)

            name = 'Target Gene List'

            yield (
                {'name': name, 'genes': ", ".join(set(chain.from_iterable(target_lists.values())))},
                next(col_strings), name, '', '', ''
            )
    except (FileNotFoundError, TypeError):
        pass
