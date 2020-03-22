import re
from itertools import chain, groupby
from operator import itemgetter, methodcaller

import pandas as pd
import pyparsing as pp
from django.db.models import Q
from django.http import JsonResponse
from django.views import View

from querytgdb.models import Analysis, AnalysisData, MetaKey

__all__ = ['OverviewView', 'OverviewAutocompleteView']

na_vals = re.compile(r'^(?:N/?A|NONE|NAN)$', flags=re.I)

name = pp.Word(pp.pyparsing_unicode.alphanums + '-_.:/\\') | pp.QuotedString('"', escChar='\\') | pp.QuotedString("'",
                                                                                                                  escChar='\\')
modname = name + pp.oneOf('= == !=') + name


def replace_none(args):
    if not args[1] or na_vals.match(args[1]):
        return (args[0], 'None', *args[2:])
    return args


def combine_name(item):
    if item[1]:
        return '{0[0]} ({0[1]})'.format(item), item[2]

    return item[0], item[2]


def format_rows(row):
    return {'id': row.pop('id'),
            'gene_id': row.pop('gene_id'),
            'gene_name': row.pop('gene_name'),
            'metadata': row}


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

        analyses = pd.DataFrame(qs.values_list('pk', 'tf__gene_id', 'tf__name').distinct().iterator(),
                                columns=['id', 'gene_id', 'gene_name'])

        # new colors in summary
        summary = {
            'tfs': list(map(combine_name,
                            analyses.groupby([
                                'gene_id',
                                'gene_name'
                            ]).count().reset_index().itertuples(index=False, name=None)))
        }

        analysisdata_qs = AnalysisData.objects.filter(analysis__in=qs).prefetch_related('key')

        analysisdata = pd.DataFrame(analysisdata_qs.values_list('analysis_id', 'key__name', 'value').iterator(),
                                    columns=['id', 'key', 'value'])

        for n, group in groupby(
                map(replace_none,
                    analysisdata.groupby([
                        'key',
                        'value'
                    ]).count().reset_index().itertuples(index=False, name=None)),
                key=itemgetter(0)):
            summary[n] = list(map(itemgetter(1, 2), group))

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
