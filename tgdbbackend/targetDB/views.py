# Create your views here.
from collections import OrderedDict
from itertools import chain

from rest_framework import viewsets
from rest_framework.decorators import detail_route
from rest_framework.response import Response

from querytgdb.models import Edges, MetaIddata, TargetDBTF
from .models import Nodes
from .serializers import EdgesValueSerializer, MetaValueSerializer, TFValueSerializer


class MetaValueDistinctViewSet(viewsets.ReadOnlyModelViewSet):
    serializer_class = MetaValueSerializer
    queryset = MetaIddata.objects.all()

    def list(self, request, *args, **kwargs):
        queryset = MetaIddata.objects.values("meta_value").distinct()
        serializer = MetaValueSerializer(queryset, many=True)
        return Response(serializer.data)

    @detail_route()
    def searchType(self, request, pk=None):
        """"""
        # uinput = request.query_params.get("uinput", '')
        pk = pk.upper()
        queryset = MetaIddata.objects.filter(meta_type=pk).values("meta_value").distinct()
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
