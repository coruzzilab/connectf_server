import math
import sys
from collections import Counter, OrderedDict
from functools import partial
from io import BytesIO
from itertools import chain, count, cycle, repeat, starmap, tee
from operator import attrgetter, itemgetter
from typing import Dict, Generator, Iterable, List, Optional, Set, Tuple, Union
from uuid import UUID

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as hierarchy
import seaborn as sns
from django.conf import settings
from django.core.cache import cache
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

from querytgdb.models import Analysis
from querytgdb.utils import clear_data, column_string, get_metadata, split_name, svg_font_adder
from querytgdb.utils.motif_enrichment.motif import AdditionalMotifData, MotifData, Region

sns.set()

MOTIFS = MotifData()
ADD_MOTIFS = AdditionalMotifData()


@MOTIFS.register
@ADD_MOTIFS.register
class Promoter2000(Region):
    name = '2000bp_promoter'
    description = '2000bp upstream promoter region'
    group = [1]


@MOTIFS.register
@ADD_MOTIFS.register
class Promoter1000(Region):
    default = True
    name = '1000bp_promoter'
    description = '1000bp upstream promoter region'
    group = [1]


@MOTIFS.register
@ADD_MOTIFS.register
class Promoter500(Region):
    name = '500bp_promoter'
    description = '500bp upstream promoter region'
    group = [1]


@MOTIFS.register
@ADD_MOTIFS.register
class FivePrimeUtr(Region):
    name = 'five_prime_UTR'
    description = "5' UTR"
    group = [2]


@MOTIFS.register
@ADD_MOTIFS.register
class Cds(Region):
    name = 'CDS'
    description = 'CDS'
    group = [3]


@MOTIFS.register
@ADD_MOTIFS.register
class Intron(Region):
    name = 'intron'
    description = 'intron'
    group = [4]


@MOTIFS.register
@ADD_MOTIFS.register
class ThreePrimeUtr(Region):
    name = 'three_prime_UTR'
    description = "3' UTR"
    group = [5]


@MOTIFS.register
@ADD_MOTIFS.register
class Exon(Region):
    name = 'exon'
    description = 'exon'
    group = [2, 3, 5]


@ADD_MOTIFS.register
class Mrna(Region):
    name = 'mrna'
    description = 'mrna'
    group = [2, 3, 4, 5]


CLUSTER_INFO = pd.read_csv(
    settings.MOTIF_CLUSTER,
    index_col=0
).fillna('').to_dict('index')


class MotifEnrichmentError(ValueError):
    pass


class NoEnrichedMotif(MotifEnrichmentError):
    pass


class MotifDict(OrderedDict):
    def __setitem__(self, key, value):
        if key in self:
            raise KeyError(f'{key} already exists')
        super().__setitem__(key, value)


def cluster_fisher(row):
    return fisher_exact((row[:2], row[2:]), alternative='greater')[1]


def get_list_enrichment(gene_list: Iterable[str],
                        motifs: Optional[List[str]],
                        annotated: pd.DataFrame,
                        ann_cluster_size: pd.DataFrame,
                        total: int) -> pd.Series:
    if motifs is None:
        list_cluster_dedup = annotated.loc[annotated.index.get_level_values(0).isin(gene_list), :]
    else:
        list_cluster_dedup = annotated.loc[
                             annotated.index.get_level_values(0).isin(gene_list) &
                             annotated.index.get_level_values(2).isin(motifs), :]
        ann_cluster_size = ann_cluster_size[ann_cluster_size.index.isin(motifs)]
    if list_cluster_dedup.empty:
        return pd.Series()

    list_cluster_size = list_cluster_dedup.groupby(level=2).sum()

    if motifs is None:
        list_cluster_sum = list_cluster_dedup.sum()
    else:
        list_cluster_sum = annotated.loc[annotated.index.get_level_values(0).isin(gene_list), :].sum()

    contingency = pd.concat([
        list_cluster_size,
        ann_cluster_size - list_cluster_size,
        list_cluster_sum - list_cluster_size,
        total - list_cluster_sum - ann_cluster_size + list_cluster_size
    ], axis=1, sort=False)

    p_values = contingency.fillna(0).apply(cluster_fisher, axis=1).sort_values()

    return p_values


def motif_enrichment(res: Dict[Tuple[str, Union[None, int]], Set[str]],
                     regions: Optional[List[str]] = None,
                     uid: Optional[Union[str, UUID]] = None,
                     motif_dict: Optional[Dict[Tuple[str, Union[None, int]], Optional[List[str]]]] = None,
                     alpha: float = 0.05,
                     show_reject: bool = True,
                     motif_data: MotifData = MOTIFS) -> Tuple[List[str], pd.DataFrame]:
    if not regions:
        regions = motif_data.default_regions
    else:
        regions = sorted(set(motif_data.regions) & set(regions) & set(motif_data.annotation.index.levels[1]),
                         key=motif_data.regions.index)

    if any(filter(lambda c: c > 1,
                  Counter(
                      chain.from_iterable(map(attrgetter('group'), map(motif_data.__getitem__, regions)))).values())):
        raise MotifEnrichmentError("Overlapping regions selected. Cannot have regions from the same group.")

    if motif_dict is None:
        motif_dict = OrderedDict(zip(res.keys(), repeat(None)))

    results: Dict[str, List[pd.Series]] = OrderedDict()

    if uid is not None:
        cached_region = cache.get_many([f'{uid}/{r}_enrich' for r in regions])
    else:
        cached_region = {}

    annotated = None

    for region in regions:
        try:
            region_enrich = cached_region[f'{uid}/{region}_enrich']
        except KeyError:
            if annotated is None:
                genes = set(chain.from_iterable(res.values()))
                annotated = motif_data.annotation
                annotated = annotated.loc[(
                                              annotated.index.get_level_values(0).isin(genes),
                                              regions,
                                              slice(None)
                                          ), :]
            region_enrich = list(map(partial(get_list_enrichment,
                                             annotated=annotated.loc[(slice(None), region, slice(None)), :],
                                             ann_cluster_size=getattr(motif_data, f'{region}_cluster_size'),
                                             total=motif_data.region_total.loc[region, :].squeeze()),
                                     res.values(),
                                     motif_dict.values()))
            if uid is not None:
                cache.set(f'{uid}/{region}_enrich', region_enrich)

        results[region] = region_enrich

    result_df = pd.concat(chain.from_iterable(zip(*results.values())), axis=1, sort=True)

    columns = list(starmap(
        lambda r, c: (*c, r),
        zip(
            cycle(regions),
            chain.from_iterable(zip(*tee(res.keys(), len(results))))
        )
    ))

    result_df.columns = columns

    # bonferroni correction
    result_df = result_df.stack()
    pvalues = multipletests(result_df, method='bonferroni')[1]
    result_df = pd.Series(pvalues, index=result_df.index)
    result_df = result_df.unstack()

    if show_reject:
        result_df = result_df[(result_df <= alpha).any(axis=1)]
    else:
        result_df[(result_df > alpha)] = np.nan
        result_df = result_df.dropna(how='all')

    if result_df.empty:
        raise NoEnrichedMotif

    return regions, result_df


def get_analysis_gene_list(uid: Union[str, UUID]) -> Dict[Tuple[str, Union[None, int]], Set[str]]:
    cached_data = cache.get_many([f'{uid}/tabular_output', f'{uid}/target_genes'])
    df = cached_data[f'{uid}/tabular_output']
    df = clear_data(df)
    res = OrderedDict((name, set(col.index[col.notnull()])) for name, col in df.iteritems())

    try:
        _, target_lists = cached_data[f'{uid}/target_genes']
        res[('Target Gene List', 0)] = set(chain.from_iterable(target_lists.values()))
    except KeyError:
        pass

    return res


def merge_cluster_info(df):
    for idx, *row in df.itertuples(name=None):
        info = {'name': idx}
        try:
            info.update(CLUSTER_INFO[idx])
        except KeyError:
            pass
        yield [info] + row


def get_column_info(res, metadata):
    columns = []

    for (tf, analysis_id), value in res.items():
        name, criterion, _uid = split_name(tf)

        data = OrderedDict([('name', name), ('filter', criterion)])
        if (tf, analysis_id) == ('Target Gene List', 0):
            data['genes'] = ",\n".join(value)
        else:
            data.update(metadata.loc[analysis_id, :].to_dict())

        columns.append(data)

    return columns


def get_motif_enrichment_json(uid: Union[str, UUID],
                              regions: Optional[List[str]] = None,
                              alpha: float = 0.05,
                              motif_data: MotifData = MOTIFS) -> Dict:
    res = get_analysis_gene_list(uid)

    analyses = Analysis.objects.prefetch_related('tf').filter(pk__in=map(itemgetter(1), res.keys()))
    metadata = get_metadata(analyses)

    regions, result_df = motif_enrichment(res, regions,
                                          uid=uid,
                                          motif_dict=None,
                                          alpha=alpha,
                                          show_reject=False,
                                          motif_data=motif_data)

    columns = get_column_info(res, metadata)

    result_df = result_df.where(pd.notnull(result_df), None)

    return {
        'columns': columns,
        'result': list(merge_cluster_info(result_df)),
        'regions': regions
    }


def get_additional_motif_enrichment_json(uid: Union[str, UUID],
                                         regions: Optional[List[str]] = None,
                                         motifs: Optional[List[Optional[List[str]]]] = None,
                                         use_default_motifs: bool = False,
                                         motif_data: MotifData = MOTIFS) -> Dict:
    res = get_analysis_gene_list(uid)

    analyses = Analysis.objects.prefetch_related('tf').filter(pk__in=map(itemgetter(1), res.keys()))
    metadata = get_metadata(analyses)

    if motifs is not None:
        if len(res) != len(motifs):
            raise MotifEnrichmentError("Should provide list of motifs for each analysis.")
    elif motifs is None and use_default_motifs:
        motifs = []
        for idx in map(itemgetter(1), res.keys()):
            try:
                motifs.append(motif_data.motifs[motif_data.motifs.str.startswith(metadata.at[idx, 'gene_id'])].tolist())
            except (KeyError, ValueError):
                motifs.append([])
    else:
        motifs = [None] * len(res)

    motif_dict = OrderedDict(zip(res.keys(), motifs))

    regions, result_df = motif_enrichment(res, regions,
                                          uid=None,
                                          motif_dict=motif_dict,
                                          alpha=1,
                                          show_reject=False,
                                          motif_data=motif_data)

    result_df.columns = pd.MultiIndex.from_tuples(result_df.columns)
    result_df = result_df.stack().unstack(level=0).dropna(how='all', axis=1)

    columns = []
    indeces = []

    for (tf, analysis_id), value in res.items():
        name, criterion, _uid = split_name(tf)

        data = OrderedDict([('name', name), ('filter', criterion)])
        motif_list = motif_dict[(tf, analysis_id)]
        data['motifs'] = motif_list

        if motif_list:
            for m in motif_list:
                indeces.append((tf, analysis_id, m))
        else:
            indeces.append((tf, analysis_id, None))

        if (tf, analysis_id) == ('Target Gene List', 0):
            data['genes'] = ",\n".join(value)
        else:
            data.update(metadata.loc[analysis_id, :].to_dict())

        columns.append(data)

    result_df = result_df.reindex(columns=pd.MultiIndex.from_tuples(indeces))
    result_df = result_df.where(pd.notnull(result_df), None)
    return {
        'columns': columns,
        'result': list(result_df.itertuples(name=None, index=True))
    }


def make_motif_enrichment_heatmap_columns(res, fields: Optional[Iterable[str]] = None) -> List[str]:
    analyses = Analysis.objects.filter(pk__in=map(itemgetter(1), res.keys())).prefetch_related('tf')

    metadata = None
    if fields:
        metadata = get_metadata(analyses, fields)

    columns = []

    for (name, analysis_id), col_name in zip(res.keys(), map(column_string, count(1))):
        try:
            tf = analyses.get(pk=analysis_id).tf
            col = '{0} — {1}{2}'.format(col_name, tf.gene_id, f' ({tf.name})' if tf.name else '')

            try:
                col += f" [{metadata.loc[analysis_id, :].str.cat(sep=', ')}]"
            except (AttributeError, KeyError):
                pass

            columns.append(col)
        except Analysis.DoesNotExist:
            columns.append('{} — {}'.format(col_name, name))

    return columns


def get_motif_enrichment_heatmap(uid: Union[str, UUID],
                                 regions: Optional[List[str]] = None,
                                 alpha: float = 0.05,
                                 lower_bound: Optional[float] = None,
                                 upper_bound: Optional[float] = None,
                                 fields: Optional[Iterable[str]] = None,
                                 motif_data: MotifData = MOTIFS) -> BytesIO:
    res = get_analysis_gene_list(uid)
    regions, result_df = motif_enrichment(res, regions, uid, alpha=alpha, motif_data=motif_data)

    result_df = result_df.clip(lower=sys.float_info.min)
    result_df = -np.log10(result_df)

    result_df = result_df.rename(
        index={idx: "{} ({})".format(idx, CLUSTER_INFO[idx]['Family']) for idx in result_df.index})

    columns = make_motif_enrichment_heatmap_columns(res, fields)

    result_df.columns = starmap(lambda r, x: f'{x} {r}',
                                zip(cycle(regions), chain.from_iterable(zip(*tee(columns, len(regions))))))

    result_df = result_df.T

    opts = {}

    rows, cols = result_df.shape

    opts['row_cluster'] = rows > 1
    opts['col_cluster'] = cols > 1

    if rows > 1:
        opts['row_linkage'] = hierarchy.linkage(result_df.values, method='average', optimal_ordering=True)
        if len(regions) > 1:
            opts['row_colors'] = [motif_data.colors[r] for r, s in zip(cycle(regions), result_df.index)]

    if cols > 1:
        opts['col_linkage'] = hierarchy.linkage(result_df.values.T, method='average', optimal_ordering=True)

    plt.figure()
    heatmap_graph = sns.clustermap(result_df,
                                   method="ward",
                                   cmap="YlGnBu",
                                   cbar_kws={'label': 'Enrichment (-log10 p)'},
                                   xticklabels=1,
                                   yticklabels=1,
                                   vmin=lower_bound,
                                   vmax=upper_bound,
                                   **opts)
    # bug in matplotlib 3.1.1
    bottom, top = heatmap_graph.ax_heatmap.get_ylim()
    heatmap_graph.ax_heatmap.set_ylim(math.ceil(bottom), math.floor(top))

    plt.setp(heatmap_graph.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    plt.setp(heatmap_graph.ax_heatmap.xaxis.get_majorticklabels(), rotation=270)
    if fields:
        heatmap_graph.ax_heatmap.set_ylabel(f'Additional fields: [{", ".join(fields)}]', rotation=270, labelpad=15)

    buffer = BytesIO()
    heatmap_graph.savefig(buffer)
    plt.close()
    buffer.seek(0)
    svg_font_adder(buffer)

    return buffer


TableRow = Tuple[Dict, str, str, str, str, str]


def get_motif_enrichment_heatmap_table(uid: Union[str, UUID]) -> Generator[TableRow, None, None]:
    cached_data = cache.get_many([f'{uid}/tabular_output', f'{uid}/target_genes'])
    df = cached_data[f'{uid}/tabular_output']
    df = clear_data(df)

    analyses = Analysis.objects.filter(pk__in=df.columns.get_level_values(1))
    metadata = get_metadata(analyses)

    col_strings = map(column_string, count(1))  # this is an infinite generator lazy evaluation only

    for (name, analysis_id), col_str in zip(df.columns, col_strings):
        name, criterion, _uid = split_name(name)

        info = OrderedDict([('name', name), ('filter', criterion)])

        try:
            try:
                gene_name = metadata.at[analysis_id, 'gene_name']
            except KeyError:
                gene_name = ''
            info.update(metadata.loc[analysis_id, :].to_dict())
        except Analysis.DoesNotExist:
            gene_name = ''

        yield (info, col_str, name, criterion, gene_name, analysis_id)

    try:
        _, target_lists = cached_data[f'{uid}/target_genes']

        name = 'Target Gene List'

        yield (
            {'name': name, 'genes': ", ".join(set(chain.from_iterable(target_lists.values())))},
            next(col_strings), name, '', '', ''
        )
    except KeyError:
        pass
