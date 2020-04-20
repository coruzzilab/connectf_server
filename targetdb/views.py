import json
import os
from itertools import filterfalse
from operator import methodcaller
from typing import Generator, Iterable, List

from django.conf import settings
from django.core.files.storage import FileSystemStorage
from django.db.models import F
from django.http import JsonResponse
from django.utils.decorators import method_decorator
from django.views import View
from django.views.decorators.csrf import ensure_csrf_cookie

from querytgdb.models import Analysis, AnalysisData, Annotation, EdgeData, EdgeType, MetaKey

gene_lists_storage = FileSystemStorage(settings.GENE_LISTS)
networks_storage = FileSystemStorage(settings.TARGET_NETWORKS)
NAMED_QUERIES = getattr(settings, 'NAMED_QUERIES', {})


def get_lists(files: Iterable) -> Generator[str, None, None]:
    for f in files:
        yield f.split(os.path.extsep, 1)[0]


def check_regulation(instance: Analysis):
    return instance.regulation_set.exists()


class TFView(View):
    @method_decorator(ensure_csrf_cookie)
    def get(self, request, *args, **kwargs):
        all_genes = request.GET.get('all', '1')

        if all_genes == '1':
            queryset = [
                *({'value': k, 'name': v} for k, v in NAMED_QUERIES.items()),
                {'value': 'multitype'},
                {'value': 'all_tfs'}
            ]
            # remove andalltfs because it is not useful beyond a handfull of datasets
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
    @method_decorator(ensure_csrf_cookie)
    def get(self, request, *args, **kwargs):
        return JsonResponse(list(EdgeType.objects.values_list("name", flat=True)), safe=False)


class InterestingListsView(View):
    @method_decorator(ensure_csrf_cookie)
    def get(self, request, *args, **kwargs):
        try:
            directories, files = gene_lists_storage.listdir('./')

            return JsonResponse(sorted(get_lists(filterfalse(methodcaller('startswith', '.'), files))), safe=False)
        except FileNotFoundError:
            return JsonResponse([], safe=False, status=404)


class InterestingNetworksView(View):
    @method_decorator(ensure_csrf_cookie)
    def get(self, request, *args, **kwargs):
        try:
            directories, files = networks_storage.listdir('./')

            return JsonResponse(sorted(get_lists(files)), safe=False)
        except FileNotFoundError:
            return JsonResponse([], safe=False)


class KeyView(View):
    @method_decorator(ensure_csrf_cookie)
    def get(self, request):
        tfs = set(filter(None, request.GET.getlist('tf')))
        try:
            all_keys = json.loads(request.GET.get('all', 'false'))
        except json.JSONDecodeError:
            all_keys = False

        if tfs & {'all_tfs', 'andalltfs', 'multitype', *NAMED_QUERIES.keys()}:
            tfs = set()

        queryset = ['targeted_by']

        if tfs:
            analyses = Analysis.objects.filter(tf__gene_id__in=tfs)

            if all_keys:
                meta_qs = MetaKey.objects.filter(
                    analysisdata__analysis__in=analyses
                )
            else:
                if EdgeData.objects.filter(tf__gene_id__in=tfs).exists():
                    queryset.append('additional_edge')
                if any(analysis.regulation_set.exists() for analysis in analyses):
                    queryset[0:0] = ['log2fc', 'pvalue']

                meta_qs = MetaKey.objects.filter(
                    analysisdata__analysis__in=analyses,
                    searchable=True
                )

            queryset.extend(
                meta_qs.distinct().values_list('name', flat=True))
        else:
            if all_keys:
                meta_qs = MetaKey.objects
            else:
                queryset[0:0] = ['log2fc', 'pvalue', 'additional_edge']
                meta_qs = MetaKey.objects.filter(searchable=True)

            queryset.extend(meta_qs.values_list('name', flat=True))

        return JsonResponse(queryset, safe=False)


class ValueView(View):
    @method_decorator(ensure_csrf_cookie)
    def get(self, request, key: str) -> JsonResponse:
        tfs = set(filter(None, request.GET.getlist('tf')))

        if tfs & {'all_tfs', 'andalltfs', 'multitype', *NAMED_QUERIES.keys()}:
            tfs = set()

        key = key.upper()

        queryset: List[str] = []

        if key in ('PVALUE', 'LOG2FC', 'TARGETED_BY'):
            return JsonResponse(queryset, safe=False)
        elif key == 'ADDITIONAL_EDGE':
            if tfs:
                queryset.extend(
                    EdgeData.objects.filter(tf__gene_id__in=tfs).distinct().values_list('type__name', flat=True))
                return JsonResponse(queryset, safe=False)

            queryset.extend(EdgeType.objects.distinct().values_list('name', flat=True))
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
