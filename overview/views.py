from itertools import chain, cycle
from operator import methodcaller

import numpy as np
import pandas as pd
import pyparsing as pp
from django.db.models import Q
from django.http import JsonResponse
from django.views import View

from querytgdb.models import Analysis, AnalysisData, MetaKey
from .utils import COLORS, TYPE_COLOR, modname, na_vals

__all__ = ['OverviewView', 'OverviewAutocompleteView']


def combine_name(item):
    if item[0][1]:
        return '{0[0]} ({0[1]})'.format(item[0]), item[1]

    return item[0][0], item[1]


def format_rows(row):
    return {'id': row.pop('id'),
            'gene_id': row.pop('gene_id'),
            'gene_name': row.pop('gene_name'),
            'metadata': row}


def get_type(s: pd.Series):
    u: np.ndarray = s.unique()

    if u.size == 1:
        return u.squeeze()

    return np.nan


def get_type_color(df: pd.DataFrame):
    if df['type'].isna().any():
        return pd.DataFrame(zip(df['count'], cycle(COLORS)), index=df.index, columns=['count', 'type'])
    return TYPE_COLOR(df)


class OverviewView(View):
    def get(self, request):
        search = request.GET.get('search')

        qs = (Analysis.objects
              .prefetch_related('tf'))

        if search:
            try:
                key, oper, value = modname.parseString(search, parseAll=True)
                q = Q(analysisdata__key__name__icontains=key)
                if oper == '!=':
                    q &= ~Q(analysisdata__value__icontains=value)
                else:
                    q &= Q(analysisdata__value__icontains=value)
                qs = qs.filter(q)
            except pp.ParseException:
                qs = qs.filter(Q(tf__gene_id__icontains=search) |
                               Q(tf__name__icontains=search) |
                               # Q(analysisdata__key__name__icontains=search) |
                               Q(analysisdata__value__icontains=search))

        analyses = pd.DataFrame(qs.values_list('pk', 'tf__gene_id', 'tf__name').distinct().iterator(chunk_size=2000),
                                columns=['id', 'gene_id', 'gene_name'])

        # new colors in summary
        summary = {
            'tfs': [(n, count, color) for (n, count), color in zip(
                map(combine_name,
                    analyses.groupby([
                        'gene_id',
                        'gene_name'
                    ]).count().itertuples(index=True, name=None)),
                cycle(COLORS))]
        }

        analysisdata_qs = AnalysisData.objects.filter(analysis__in=qs).prefetch_related('key')

        analysisdata = pd.DataFrame(analysisdata_qs.values_list('analysis_id', 'key__name', 'value').iterator(chunk_size=2000),
                                    columns=['id', 'key', 'value'])

        analysisdata['value'] = analysisdata['value'].str.replace(na_vals, 'None', regex=True)

        exp_type = analysisdata.loc[analysisdata['key'] == 'EXPERIMENT_TYPE', ['id', 'value']]
        exp_type['value'] = exp_type['value'].str.upper()

        analysisdata_group = analysisdata.merge(exp_type, on='id')
        analysisdata_group = (analysisdata_group
                              .groupby(['key', 'value_x'])
                              .agg({'value_x': 'count', 'value_y': get_type})
                              .rename(columns={'value_x': 'count', 'value_y': 'type'}))

        analysisdata_group = analysisdata_group.groupby(level=0).apply(get_type_color)

        for n, group in analysisdata_group.groupby(level=0):
            summary[n] = list(group.reset_index(level=1).itertuples(index=False, name=None))

        return JsonResponse({
            'datasets': list(map(format_rows,
                                 (analyses
                                  .merge(
                                     analysisdata.pivot(index='id', columns='key', values='value').sort_index(axis=1),
                                     left_on='id', right_index=True, how='left')
                                  .fillna('')
                                  .sort_values(['gene_id', 'id'])
                                  .to_dict(orient='records')))),
            'summary': summary,
            'total': Analysis.objects.count()
        })


class OverviewAutocompleteView(View):
    def get(self, request):
        analyses = Analysis.objects.prefetch_related('tf', 'analysisdata_set').distinct()

        values = sorted(
            filter(None, chain(
                *zip(*analyses.values_list('tf__name', 'tf__gene_id')),
                analyses.values_list('analysisdata__value', flat=True),
                MetaKey.objects.distinct().values_list('name', flat=True)
            )),
            key=methodcaller('lower'))

        return JsonResponse(values, safe=False)
