import os
from typing import Generator, Iterable, List

from django.conf import settings
from django.core.files.storage import FileSystemStorage
from django.db.models import F
from django.http import JsonResponse
from django.views import View

from querytgdb.models import Analysis, AnalysisData, Annotation, EdgeData, EdgeType, MetaKey

gene_lists_storage = FileSystemStorage(settings.GENE_LISTS)
networks_storage = FileSystemStorage(settings.TARGET_NETWORKS)


def get_lists(files: Iterable) -> Generator[str, None, None]:
    for f in files:
        yield f.split(os.path.extsep, 1)[0]


def check_regulation(instance: Analysis):
    return instance.regulation_set.exists()


class TFView(View):
    def get(self, request, *args, **kwargs):
        all_genes = request.GET.get('all', '1')

        if all_genes == '1':
            queryset = [{'value': 'oralltfs'},
                        # remove andalltfs because it is not useful beyond a handfull of datasets
                        # OrderedDict([('gene_id', 'andalltfs')]),
                        {'value': 'multitype'}]
        else:
            queryset = []

        queryset.extend(
            Annotation.objects.filter(analysis__in=Analysis.objects.all())
                .distinct()
                .annotate(value=F('gene_id'))
                .values('value', 'name')
                .order_by('value', 'name')
                .iterator())

        return JsonResponse(queryset, safe=False)


class EdgeListView(View):
    def get(self, request, *args, **kwargs):
        return JsonResponse(list(EdgeType.objects.values_list("name", flat=True)), safe=False)


class InterestingListsView(View):
    def get(self, request, *args, **kwargs):
        directories, files = gene_lists_storage.listdir('./')

        return JsonResponse(list(get_lists(files)), safe=False)


class InterestingNetworksView(View):
    def get(self, request, *args, **kwargs):
        directories, files = networks_storage.listdir('./')

        return JsonResponse(list(get_lists(files)), safe=False)


class KeyView(View):
    def get(self, request):
        tfs = set(request.GET.getlist('tf'))

        if tfs & {'oralltfs', 'andalltfs', 'multitype'}:
            tfs = set()

        queryset = ['has_column']

        if tfs:
            if EdgeData.objects.filter(tf__gene_id__in=tfs).exists():
                queryset.append('additional_edge')

            analyses = Analysis.objects.filter(tf__gene_id__in=tfs)

            if any(analysis.regulation_set.exists() for analysis in analyses):
                queryset[0:0] = ['fc', 'pvalue']

            queryset.extend(
                MetaKey.objects.filter(
                    analysisdata__analysis__in=analyses,
                    searchable=True
                ).distinct().values_list('name', flat=True))
        else:
            queryset[0:0] = ['fc', 'pvalue', 'additional_edge']
            queryset.extend(MetaKey.objects.filter(searchable=True).values_list('name', flat=True))

        return JsonResponse(queryset, safe=False)


class ValueView(View):
    def get(self, request, key: str) -> JsonResponse:
        tfs = set(request.GET.getlist('tf'))

        if tfs & {'oralltfs', 'andalltfs', 'multitype'}:
            tfs = set()

        key = key.upper()

        queryset: List[str] = []

        if key in ('PVALUE', 'FC'):
            return JsonResponse(queryset, safe=False)
        elif key == 'ADDITIONAL_EDGE':
            if tfs:
                queryset.extend(
                    EdgeData.objects.filter(tf__gene_id__in=tfs).distinct().values_list('type__name', flat=True))
                return JsonResponse(queryset, safe=False)

            queryset.extend(EdgeType.objects.distinct().values_list('name', flat=True))
            return JsonResponse(queryset, safe=False)
        elif key == 'HAS_COLUMN':
            queryset.append('EDGE')

            if tfs:
                if any(map(check_regulation, Analysis.objects.filter(tf__gene_id__in=tfs))):
                    queryset.extend(('Pvalue', 'FC'))
                if EdgeData.objects.filter(tf__gene_id__in=tfs).exists():
                    queryset.append('additional_edge')
            else:
                queryset.extend(('Pvalue', 'FC', 'additional_edge'))

            return JsonResponse(queryset, safe=False)
        else:
            if tfs:
                analyses = Analysis.objects.filter(tf__gene_id__in=tfs)

                queryset.extend(AnalysisData.objects.filter(
                    analysis__in=analyses,
                    key__name__iexact=key
                ).distinct().values_list('value', flat=True))
            else:
                queryset.extend(
                    AnalysisData.objects.filter(key__name__iexact=key).distinct().values_list('value', flat=True))

            return JsonResponse(queryset, safe=False)
