# Create your views here.
from collections import OrderedDict
from itertools import chain

from rest_framework import viewsets
from rest_framework.decorators import detail_route
from rest_framework.response import Response

from querytgdb.models import AnalysisIddata, Edges, MetaIddata, Metadata, TargetDBTF
from .models import Nodes
from .serializers import EdgesValueSerializer, ExperimentIdSerializer, MetaValueSerializer, TFValueSerializer


class MetaValueDistinctViewSet(viewsets.ReadOnlyModelViewSet):
    serializer_class = ExperimentIdSerializer
    queryset = Metadata.objects.values("meta_fullid").distinct()

    @detail_route()
    def searchType(self, request, pk=None):
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
    serializer_class = TFValueSerializer
    queryset = Nodes.objects.all().values("text").distinct()

    @detail_route()
    def searchName(self, request, pk=None):
        # uinput = request.query_params.get("uinput", None)
        queryset = TargetDBTF.objects.raw("SELECT db_tf_id, db_tf_agi, ath_name FROM querytgdb_targetdbtf "
                                          "LEFT JOIN querytgdb_annotation ON agi_id = db_tf_agi")

        serializer = TFValueSerializer(
            chain([OrderedDict([('db_tf_agi', 'OR [ALLTF]')]), OrderedDict([('db_tf_agi', 'AND [ALLTF]')])], queryset),
            many=True)

        return Response(serializer.data)


class EdgesValueDistinctViewSet(viewsets.ReadOnlyModelViewSet):
    serializer_class = EdgesValueSerializer
    queryset = Edges.objects.values("edge_name").distinct()

    @detail_route()
    def searchName(self, request, pk=None):
        queryset = Edges.objects.values("edge_name").distinct()
        serializer = EdgesValueSerializer(queryset, many=True)
        return Response(serializer.data)
