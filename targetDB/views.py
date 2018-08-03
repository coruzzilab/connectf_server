import os
from collections import OrderedDict
from typing import Generator, Iterable

from django.core.files.storage import FileSystemStorage
from django.db import connection
from rest_framework import views
from rest_framework.response import Response

from querytgdb.models import AnalysisIddata, EdgeData, EdgeType, Edges, Interactions, MetaIddata, ReferenceId
from .serializers import TFValueSerializer

storage = FileSystemStorage('commongenelists/')


def get_lists(files: Iterable) -> Generator[str, None, None]:
    for f in files:
        name, ext = os.path.splitext(f)
        if ext == '.txt':
            yield name


def check_regulation(instance: ReferenceId):
    return instance.regulation_set.exists()


class TFView(views.APIView):
    def get(self, request, *args, **kwargs):
        queryset = [OrderedDict([('db_tf_agi', 'oralltfs')]),
                    OrderedDict([('db_tf_agi', 'andalltfs')])]

        with connection.cursor() as cursor:
            cursor.execute("SELECT DISTINCT db_tf_id, db_tf_agi, ath_name FROM querytgdb_targetdbtf "
                           "LEFT JOIN querytgdb_annotation ON agi_id = db_tf_agi "
                           "ORDER BY db_tf_agi, ath_name")

            for db_tf_id, db_tf_agi, ath_name in cursor:
                queryset.append(OrderedDict([
                    ("db_tf_id", db_tf_id),
                    ("db_tf_agi", db_tf_agi),
                    ("ath_name", ath_name)
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

        if tfs & {'oralltfs', 'andalltfs'}:
            tfs = set()

        queryset = ['pvalue', 'edge', 'fc', 'has_column']

        if tfs:
            refs = Interactions.objects.filter(
                db_tf_id__db_tf_agi__in=tfs).distinct().values_list('ref_id_id', flat=True)

            if EdgeData.objects.filter(tf__agi_id__in=tfs).exists():
                queryset.append('edge_properties')

            queryset.extend(AnalysisIddata.objects.filter(
                analysis_id__referenceid__ref_id__in=refs
            ).distinct().values_list('analysis_type', flat=True))

            queryset.extend(MetaIddata.objects.filter(
                meta_id__referenceid__ref_id__in=refs
            ).distinct().values_list('meta_type', flat=True))
        else:
            queryset.append('edge_properties')
            queryset.extend(AnalysisIddata.objects.distinct().values_list('analysis_type', flat=True))
            queryset.extend(MetaIddata.objects.distinct().values_list('meta_type', flat=True))

        return Response(queryset)


class ValueView(views.APIView):
    def get(self, request, key: str) -> Response:
        tfs = set(request.GET.getlist('tf'))

        if tfs & {'oralltfs', 'andalltfs'}:
            tfs = set()

        key = key.upper()

        if key in ('PVALUE', 'FC'):
            return Response([])
        elif key == 'EDGE':
            if tfs:
                return Response(Interactions.objects.filter(db_tf_id__db_tf_agi__in=tfs).distinct().values_list(
                    'edge_id__edge_name', flat=True))
            return Response(Edges.objects.distinct().values_list('edge_name', flat=True))
        elif key == 'EDGE_PROPERTIES':
            if tfs:
                return Response(
                    EdgeData.objects.filter(tf__agi_id__in=tfs).distinct().values_list('type__name', flat=True))
            return Response(EdgeType.objects.distinct().values_list('name', flat=True))
        elif key == 'HAS_COLUMN':
            queryset = ['EDGE']

            if tfs:
                if any(map(check_regulation,
                           ReferenceId.objects.filter(interactions__db_tf_id__db_tf_agi__in=tfs).distinct())):
                    queryset.extend(('Pvalue', 'FC'))
                if EdgeData.objects.filter(tf__agi_id__in=tfs).exists():
                    queryset.append('edge_properties')
            else:
                queryset.extend(('Pvalue', 'FC', 'edge_properties'))

            return Response(queryset)
        else:
            queryset = []

            if tfs:
                refs = Interactions.objects.filter(
                    db_tf_id__db_tf_agi__in=tfs).distinct().values_list('ref_id_id', flat=True)

                queryset.extend(AnalysisIddata.objects.filter(
                    analysis_id__referenceid__ref_id__in=refs,
                    analysis_type__iexact=key).distinct().values_list(
                    'analysis_value',
                    flat=True))

                queryset.extend(
                    MetaIddata.objects.filter(
                        meta_id__referenceid__ref_id__in=refs,
                        meta_type__iexact=key).distinct().values_list('meta_value', flat=True))
            else:
                queryset.extend(
                    AnalysisIddata.objects.filter(analysis_type__iexact=key).distinct().values_list('analysis_value',
                                                                                                    flat=True))
                queryset.extend(
                    MetaIddata.objects.filter(meta_type__iexact=key).distinct().values_list('meta_value', flat=True))

            return Response(queryset)
