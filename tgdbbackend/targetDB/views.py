# Create your views here.

from collections import OrderedDict
from itertools import chain

import pandas as pd
from django.db import connection
from rest_framework import viewsets
from rest_framework.decorators import detail_route, list_route
from rest_framework.response import Response

from querytgdb.models import AnalysisIddata, Annotation, Edges, MetaIddata, Metadata
from .serializers import AnnotationSerializer, EdgesValueSerializer, ExperimentIdSerializer, MetaValueSerializer, \
    TFValueSerializer


class MetaValueDistinctViewSet(viewsets.ReadOnlyModelViewSet):
    serializer_class = ExperimentIdSerializer
    queryset = Metadata.objects.values("meta_fullid").distinct()

    @detail_route()
    def search_type(self, request, pk=None):
        """"""
        # uinput = request.query_params.get("uinput", '')
        pk = pk.upper()
        queryset = chain(
            MetaIddata.objects.filter(meta_type=pk).values("meta_value").distinct(),
            AnalysisIddata.objects.filter(analysis_type=pk).extra(
                select={'meta_value': 'analysis_value'}).values('meta_value').distinct()
        )
        serializer = MetaValueSerializer(queryset, many=True)
        return Response(serializer.data)


class TFValueDistinctViewSet(viewsets.ReadOnlyModelViewSet):
    serializer_class = AnnotationSerializer
    queryset = Annotation.objects.filter(ath_gene_type='TXNFACTOR')

    @list_route()
    def search_name(self, request, *args, **kwargs):
        all_opt = request.GET.get('all')

        queryset = []

        with connection.cursor() as cursor:
            cursor.execute("SELECT DISTINCT db_tf_id, db_tf_agi, ath_name, meta_fullid FROM querytgdb_targetdbtf "
                           "LEFT JOIN querytgdb_annotation ON agi_id = db_tf_agi "
                           "INNER JOIN querytgdb_interactions ON db_tf_id = db_tf_id_id "
                           "INNER JOIN querytgdb_referenceid ON ref_id_id = ref_id "
                           "INNER JOIN querytgdb_metadata ON meta_id_id = meta_id "
                           "ORDER BY db_tf_agi, ath_name, meta_fullid")
            data = pd.DataFrame.from_records(
                list(cursor),
                columns=("db_tf_id", "db_tf_agi", "ath_name", "meta_fullid"))

            for (db_tf_id, db_tf_agi, ath_name), group in data.groupby(["db_tf_id", "db_tf_agi", "ath_name"]):
                queryset.append(OrderedDict([
                    ("db_tf_id", db_tf_id),
                    ("db_tf_agi", db_tf_agi),
                    ("ath_name", ath_name),
                    ("meta_fullid", group.meta_fullid.tolist())
                ]))

        if all_opt:
            queryset = [OrderedDict([('db_tf_agi', 'OR [ALLTF]')]),
                        OrderedDict([('db_tf_agi', 'AND [ALLTF]')])] + queryset

        serializer = TFValueSerializer(queryset, many=True)

        return Response(serializer.data)


class EdgesValueDistinctViewSet(viewsets.ReadOnlyModelViewSet):
    serializer_class = EdgesValueSerializer
    queryset = Edges.objects.values("edge_name").distinct()

    @detail_route()
    def search_name(self, request, *args, **kwargs):
        queryset = Edges.objects.values("edge_name").distinct()
        serializer = EdgesValueSerializer(queryset, many=True)
        return Response(serializer.data)
