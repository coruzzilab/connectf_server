import math
import sys
from collections import Counter, OrderedDict
from functools import partial
from io import BytesIO
from itertools import chain, count, cycle, starmap, tee
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
from querytgdb.utils import annotations, clear_data, column_string, get_metadata, split_name, svg_font_adder
from querytgdb.utils.motif_enrichment.motif import MotifData, Region

sns.set()

MOTIF = MotifData(settings.MOTIF_ANNOTATION)


@MOTIF.register
class Promoter2000(Region):
    name = '2000bp_promoter'
    description = '2000bp upstream promoter region'
    group = [1]


@MOTIF.register
class Promoter1000(Region):
    default = True
    name = '1000bp_promoter'
    description = '1000bp upstream promoter region'
    group = [1]


@MOTIF.register
class Promoter500(Region):
    name = '500bp_promoter'
    description = '500bp upstream promoter region'
    group = [1]


@MOTIF.register
class FivePrimeUtr(Region):
    name = 'five_prime_UTR'
    description = "5' UTR"
    group = [2]


@MOTIF.register
class Cds(Region):
    name = 'CDS'
    description = 'CDS'
    group = [3]


@MOTIF.register
class Intron(Region):
    name = 'intron'
    description = 'intron'
    group = [4]


@MOTIF.register
class ThreePrimeUtr(Region):
    name = 'three_prime_UTR'
    description = "3' UTR"
    group = [5]


@MOTIF.register
class Exon(Region):
    name = 'exon'
    description = 'exon'
    group = [2, 3, 5]


CLUSTER_INFO = pd.read_csv(
    settings.MOTIF_CLUSTER,
    index_col=0
).fillna('').to_dict('index')

COLORS = dict(zip(MOTIF.regions, sns.color_palette("husl", len(MOTIF.regions))))


class MotifEnrichmentError(ValueError):
    pass


class NoEnrichedMotif(MotifEnrichmentError):
    pass


def cluster_fisher(row):
    return fisher_exact((row[:2], row[2:]), alternative='greater')[1]


def get_list_enrichment(gene_list: Iterable[str],
                        annotated: pd.DataFrame,
                        ann_cluster_size: pd.DataFrame) -> pd.Series:
    list_cluster_dedup = annotated.loc[annotated.index.get_level_values(0).isin(gene_list), :]
    list_cluster_size = list_cluster_dedup.groupby(level=2).sum()
    list_cluster_sum = list_cluster_dedup.sum()

    p_values = pd.concat([
        list_cluster_size,
        ann_cluster_size - list_cluster_size,
        list_cluster_sum - list_cluster_size,
        annotated.sum() - list_cluster_sum - ann_cluster_size + list_cluster_size
    ], axis=1, sort=False).fillna(0).apply(cluster_fisher, axis=1).sort_values()

    return p_values


def motif_enrichment(res: Dict[Tuple[str, Union[None, int]], Set[str]],
                     regions: Optional[List[str]] = None,
                     uid: Optional[Union[str, UUID]] = None,
                     alpha: float = 0.05,
                     show_reject: bool = True) -> Tuple[List[str], pd.DataFrame]:
    if not regions:
        regions = MOTIF.default_regions
    else:
        regions = sorted(set(MOTIF.regions) & set(regions), key=MOTIF.regions.index)

    if any(filter(lambda c: c > 1,
                  Counter(chain.from_iterable(map(attrgetter('group'), map(MOTIF.__getitem__, regions)))).values())):
        raise MotifEnrichmentError("Overlapping regions selected. Cannot have regions from the same group.")

    results: Dict[str, List[pd.Series]] = OrderedDict()

    cached_region = cache.get_many([f'{uid}/{r}_enrich' for r in regions])

    for region in regions:
        try:
            region_enrich = cached_region[f'{uid}/{region}_enrich']
        except KeyError:
            region_enrich = list(map(partial(get_list_enrichment,
                                             annotated=getattr(MOTIF, region),
                                             ann_cluster_size=getattr(MOTIF, f'{region}_cluster_size')),
                                     res.values()))
            cache.set(f'{uid}/{region}_enrich', region_enrich)

        results[region] = region_enrich

    df = pd.concat(chain.from_iterable(zip(*results.values())), axis=1, sort=True)

    columns = list(starmap(
        lambda r, c: c + (r,),
        zip(
            cycle(regions),
            chain.from_iterable(zip(*tee(res.keys(), len(results))))
        )
    ))

    df.columns = columns

    # bonferroni correction
    df = df.stack()
    pvalues = multipletests(df, method='bonferroni')[1]
    df = pd.Series(pvalues, index=df.index)
    df = df.unstack()

    if show_reject:
        df = df[(df <= alpha).any(axis=1)]
    else:
        df[df > alpha] = np.nan
        df = df.dropna(how='all')

    if df.empty:
        raise NoEnrichedMotif

    return regions, df


def get_analysis_gene_list(uid: Union[str, UUID]) -> Dict[Tuple[str, Union[None, int]], Set[str]]:
    cached_data = cache.get_many([f'{uid}/tabular_output', f'{uid}/target_genes'])
    df = cached_data[f'{uid}/tabular_output']
    df = clear_data(df)
    res = OrderedDict((name, set(col.index[col.notnull()])) for name, col in df.iteritems())

    try:
        _, target_lists = cached_data[f'{uid}/target_genes']
        res[('Target Gene List', None)] = set(chain.from_iterable(target_lists.values()))
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


def get_motif_enrichment_json(uid: Union[str, UUID],
                              regions: Optional[List[str]] = None,
                              alpha: float = 0.05) -> Dict:
    res = get_analysis_gene_list(uid)
    regions, df = motif_enrichment(res, regions, uid, alpha=alpha, show_reject=False)

    analyses = Analysis.objects.prefetch_related('tf').filter(pk__in=map(itemgetter(1), res.keys()))
    metadata = get_metadata(analyses)

    columns = []

    for (tf, analysis_id), value in res.items():
        name, criterion, _uid = split_name(tf)

        data = OrderedDict([('name', name), ('filter', criterion)])
        try:
            data.update(metadata.loc[analysis_id, :].to_dict())
        except KeyError:
            if (tf, analysis_id) == ('Target Gene List', None):
                data['genes'] = ",\n".join(value)
            else:
                raise

        columns.append(data)

    df = df.where(pd.notnull(df), None)

    return {
        'columns': columns,
        'result': list(merge_cluster_info(df)),
        'regions': regions
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
                                 fields: Optional[Iterable[str]] = None) -> BytesIO:
    res = get_analysis_gene_list(uid)
    regions, df = motif_enrichment(res, regions, uid, alpha=alpha)

    df = df.clip(lower=sys.float_info.min)
    df = -np.log10(df)

    df = df.rename(index={idx: "{} ({})".format(idx, CLUSTER_INFO[idx]['Family']) for idx in df.index})

    columns = make_motif_enrichment_heatmap_columns(res, fields)

    df.columns = starmap(lambda r, x: f'{x} {r}',
                         zip(cycle(regions), chain.from_iterable(zip(*tee(columns, len(regions))))))

    df = df.T

    opts = {}

    rows, cols = df.shape

    opts['row_cluster'] = rows > 1
    opts['col_cluster'] = cols > 1

    if rows > 1:
        opts['row_linkage'] = hierarchy.linkage(df.values, method='average', optimal_ordering=True)
        if len(regions) > 1:
            opts['row_colors'] = [COLORS[r] for r, s in zip(cycle(regions), df.index)]

    if cols > 1:
        opts['col_linkage'] = hierarchy.linkage(df.values.T, method='average', optimal_ordering=True)

    plt.figure()
    heatmap_graph = sns.clustermap(df,
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

    col_strings = map(column_string, count(1))  # this is an infinite generator lazy evaluation only

    for (name, analysis_id), col_str in zip(df.columns, col_strings):
        name, criterion, _uid = split_name(name)

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
        _, target_lists = cached_data[f'{uid}/target_genes']

        name = 'Target Gene List'

        yield (
            {'name': name, 'genes': ", ".join(set(chain.from_iterable(target_lists.values())))},
            next(col_strings), name, '', '', ''
        )
    except KeyError:
        pass
