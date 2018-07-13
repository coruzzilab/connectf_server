import os
from collections import OrderedDict
from typing import Generator, Iterable

from django.core.files.storage import FileSystemStorage
from django.db import connection
from rest_framework import views
from rest_framework.response import Response

from querytgdb.models import AnalysisIddata, Edges, MetaIddata
from .serializers import TFValueSerializer

storage = FileSystemStorage('commongenelists/')


def get_lists(files: Iterable) -> Generator[str, None, None]:
    for f in files:
        name, ext = os.path.splitext(f)
        if ext == '.txt':
            yield name


class TFView(views.APIView):
    def get(self, request, *args, **kwargs):
        queryset = [OrderedDict([('db_tf_agi', 'oralltf')]),
                    OrderedDict([('db_tf_agi', 'andalltf')])]

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


class InterestingListsView(views.APIView):
    def get(self, request, *args, **kwargs):
        directories, files = storage.listdir('./')

        return Response(get_lists(files))


class KeyView(views.APIView):
    def get(self, request):
        queryset = ['pvalue', 'edge', 'fc']
        queryset.extend(AnalysisIddata.objects.distinct().values_list('analysis_type', flat=True))
        queryset.extend(MetaIddata.objects.distinct().values_list('meta_type', flat=True))

        return Response(queryset)


class ValueView(views.APIView):
    def get(self, request, key: str) -> Response:
        key = key.upper()

        if key in ('PVALUE', 'FC'):
            return Response([])
        elif key == 'EDGE':
            return Response(Edges.objects.distinct().values_list('edge_name', flat=True))
        else:
            queryset = []

            queryset.extend(
                AnalysisIddata.objects.filter(analysis_type__iexact=key).distinct().values_list('analysis_value',
                                                                                                flat=True))
            queryset.extend(
                MetaIddata.objects.filter(meta_type__iexact=key).distinct().values_list('meta_value', flat=True))

            return Response(queryset)
