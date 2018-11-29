import os
from collections import OrderedDict
from typing import Generator, Iterable

from django.conf import settings
from django.core.files.storage import FileSystemStorage
from django.db import connection
from rest_framework import views
from rest_framework.response import Response

from querytgdb.models import Analysis, AnalysisData, EdgeData, EdgeType, MetaKey
from .serializers import TFValueSerializer

storage = FileSystemStorage(settings.GENE_LISTS)


def get_lists(files: Iterable) -> Generator[str, None, None]:
    for f in files:
        name, ext = os.path.splitext(f)
        if ext == '.txt':
            yield name


def check_regulation(instance: Analysis):
    return instance.regulation_set.exists()


class TFView(views.APIView):
    def get(self, request, *args, **kwargs):
        all_genes = request.GET.get('all', '1')

        if all_genes == '1':
            queryset = [OrderedDict([('gene_id', 'oralltfs')]),
                        # remove andalltfs because it is not use beyond a handfull of datasets
                        # OrderedDict([('gene_id', 'andalltfs')]),
                        OrderedDict([('gene_id', 'multitype')])]
        else:
            queryset = []

        with connection.cursor() as cursor:
            cursor.execute(
                "SELECT DISTINCT a.gene_id as gene_id, a.name as gene_name FROM querytgdb_analysis as e "
                "LEFT JOIN querytgdb_annotation as a ON e.tf_id = a.id "
                "ORDER BY gene_id, gene_name")

            for gene_id, gene_name in cursor:
                queryset.append(OrderedDict([
                    ("gene_id", gene_id),
                    ("gene_name", gene_name)
                ]))

        serializer = TFValueSerializer(queryset, many=True)

        return Response(serializer.data)


class EdgeListView(views.APIView):
    def get(self, request, *args, **kwargs):
        return Response(EdgeType.objects.values_list("name", flat=True))


class InterestingListsView(views.APIView):
    def get(self, request, *args, **kwargs):
        directories, files = storage.listdir('./')

        return Response(get_lists(files))


class KeyView(views.APIView):
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

        return Response(queryset)


class ValueView(views.APIView):
    def get(self, request, key: str) -> Response:
        tfs = set(request.GET.getlist('tf'))

        if tfs & {'oralltfs', 'andalltfs', 'multitype'}:
            tfs = set()

        key = key.upper()

        if key in ('PVALUE', 'FC'):
            return Response([])
        elif key == 'ADDITIONAL_EDGE':
            if tfs:
                return Response(
                    EdgeData.objects.filter(tf__gene_id__in=tfs).distinct().values_list('type__name', flat=True))
            return Response(EdgeType.objects.distinct().values_list('name', flat=True))
        elif key == 'HAS_COLUMN':
            queryset = ['EDGE']

            if tfs:
                if any(map(check_regulation, Analysis.objects.filter(tf__gene_id__in=tfs))):
                    queryset.extend(('Pvalue', 'FC'))
                if EdgeData.objects.filter(tf__gene_id__in=tfs).exists():
                    queryset.append('additional_edge')
            else:
                queryset.extend(('Pvalue', 'FC', 'additional_edge'))

            return Response(queryset)
        else:
            queryset = []

            if tfs:
                analyses = Analysis.objects.filter(tf__gene_id__in=tfs)

                queryset.extend(AnalysisData.objects.filter(
                    analysis__in=analyses,
                    key__name__iexact=key
                ).distinct().values_list('value', flat=True))
            else:
                queryset.extend(
                    AnalysisData.objects.filter(key__name__iexact=key).distinct().values_list('value', flat=True))

            return Response(queryset)
