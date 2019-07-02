import io
import sys
from itertools import count, product
from operator import itemgetter
from typing import Any, Dict, List, Optional, Union
from uuid import UUID

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as hierarchy
import seaborn as sns
from django.core.cache import cache
from scipy.stats import fisher_exact

from querytgdb.utils import annotations
from ..models import Analysis, AnalysisData
from ..utils import clear_data, column_string, split_name, get_metadata

sns.set()


def scale_df(df: pd.DataFrame) -> pd.DataFrame:
    df = df.clip(lower=sys.float_info.min)
    return -np.log10(df)


def draw_heatmap(df: pd.DataFrame, **kwargs):
    row_num, col_num = df.shape
    opts = {}
    if row_num > 1:
        opts['row_linkage'] = hierarchy.linkage(df.values, method='average', optimal_ordering=True)
    if col_num > 1:
        opts['col_linkage'] = hierarchy.linkage(df.values.T, method='average', optimal_ordering=True)

    sns_heatmap = sns.clustermap(df,
                                 cmap="YlGnBu",
                                 cbar_kws={'label': 'Enrichment (-log10 p)'},
                                 xticklabels=1,
                                 yticklabels=1,
                                 row_cluster=row_num > 1,
                                 col_cluster=col_num > 1,
                                 **kwargs,
                                 **opts)
    plt.setp(sns_heatmap.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    # plt.setp(sns_heatmap.ax_heatmap.yaxis.get_label(), text="Analyses", rotation=270)
    plt.setp(sns_heatmap.ax_heatmap.xaxis.get_majorticklabels(), rotation=270)
    plt.setp(sns_heatmap.ax_heatmap.xaxis.get_label(), text="Target Genes")

    return sns_heatmap


def gene_list_enrichment(uid: Union[str, UUID], background: Optional[int] = None, draw=True, legend=False,
                         upper=None, lower=None, fields: Optional[List[str]] = None):
    """
    Draws gene list enrichment heatmap from cached query

    :param uid:
    :param background:
    :param draw:
    :param legend:
    :param upper:
    :param lower:
    :param fields:
    :return:
    """
    # raising exception here if target genes are not uploaded by the user
    cached_data = cache.get_many([f'{uid}/target_genes', f'{uid}/tabular_output_unfiltered'])

    try:
        name_to_list, list_to_name = cached_data[f'{uid}/target_genes']
    except KeyError as e:
        raise ValueError('No target genes uploaded') from e

    try:
        query_result = cached_data[f'{uid}/tabular_output_unfiltered']
    except KeyError as e:
        raise ValueError('Query result unavailable') from e

    if background is None:
        if not annotations().empty:
            background = annotations().shape[0]
        else:
            background = 28775

    # default background:
    # A. thaliana columbia tair10 genome(28775 genes(does not include transposable elements and pseudogenes))

    query_result = clear_data(query_result)

    analyses = Analysis.objects.filter(
        pk__in=query_result.columns.get_level_values(1)
    ).distinct().prefetch_related('tf')

    # save targets for each TF into a dict
    targets = {}

    for (name, analysis_id), column in query_result.iteritems():
        analysis_targets = set(column.dropna().index.str.upper())

        name, criterion = split_name(name)

        targets[(name, criterion, analysis_id, len(analysis_targets))] = analysis_targets

    list_enrichment_pvals = pd.DataFrame(index=targets.keys(), columns=list_to_name.keys(), dtype=np.float64)

    if not legend:
        list_enrichment_count = pd.DataFrame(index=targets.keys(), columns=list_to_name.keys(), dtype=np.int_)

        colnames = {}

        for (analysis_name, analysis_list), (name, user_list) in product(targets.items(), list_to_name.items()):
            intersect_len = len(analysis_list & user_list)
            colnames[name] = "{} ({})".format(name, len(user_list))

            odds, pvalue = fisher_exact([
                [intersect_len, len(user_list - analysis_list)],
                [len(analysis_list - user_list), background - len(user_list | analysis_list)]
            ], alternative='greater')
            list_enrichment_pvals.at[analysis_name, name] = pvalue
            list_enrichment_count.at[analysis_name, name] = intersect_len

        list_enrichment_pvals = list_enrichment_pvals.rename(columns=colnames)
        list_enrichment_count = list_enrichment_count.rename(columns=colnames)

        if not draw:
            return list_enrichment_pvals, list_enrichment_count

    if draw or legend:
        orig_index = list(zip(list_enrichment_pvals.index, map(column_string, count(1))))

        metadata = None

        if fields:
            metadata = get_metadata(analyses, fields)

        indices = []

        for (name, criterion, analysis_id, l), col_name in orig_index:
            tf = analyses.get(pk=analysis_id).tf
            idx = '{} — {}{} ({})'.format(col_name, tf.gene_id, f' {tf.name}' if tf.name else '', l)

            try:
                idx += f" [{metadata.loc[analysis_id, :].fillna('None').str.cat(sep=', ')}]"
            except (KeyError, AttributeError):
                pass

            indices.append(idx)

        list_enrichment_pvals.index = indices

        if legend:
            result = []

            for idx, col_label in orig_index:
                name, criterion, analysis_id, l = idx

                info = {'name': name, 'filter': criterion}
                analysis = analyses.get(pk=analysis_id)

                gene_id = analysis.tf.gene_id

                try:
                    gene_name = annotations().at[gene_id, 'Name']
                except KeyError:
                    gene_name = ''

                info.update(analysis.analysisdata_set.values_list('key__name', 'value').iterator())

                result.append((info, col_label, name, criterion, l, gene_name, analysis_id))

            return result
        else:
            scaled_pvals = scale_df(list_enrichment_pvals)
            sns_heatmap = draw_heatmap(scaled_pvals, vmin=lower, vmax=upper)

            if fields:
                sns_heatmap.ax_heatmap.set_ylabel(f'Addational fields: [{", ".join(fields)}]',
                                                  rotation=270,
                                                  labelpad=15)

            buff = io.BytesIO()

            sns_heatmap.savefig(buff)
            plt.close()

            buff.seek(0)

            return buff


def gene_list_enrichment_json(uid) -> Dict[str, Any]:
    data, counts = gene_list_enrichment(
        uid,
        draw=False,
        legend=False
    )
    analyses = Analysis.objects.filter(pk__in=map(itemgetter(2), data.index)).prefetch_related('analysisdata_set',
                                                                                               'analysisdata_set__key',
                                                                                               'tf')

    result = []

    for ((name, criterion, analysis_id, l), *row), count_row in zip(data.itertuples(name=None),
                                                                    counts.itertuples(name=None, index=False)):
        info = {'filter': criterion, 'targets': l}
        try:
            analysis = analyses.get(pk=analysis_id)
            info['name'] = analysis.tf.gene_name_symbol
            info.update(analysis.meta_dict)
        except Analysis.DoesNotExist:
            pass

        result.append({'info': info,
                       'p-value': row,
                       'count': count_row})

    return {
        'columns': data.columns,
        'result': result
    }
