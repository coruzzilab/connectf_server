# Create your views here.
from django.db.models import Q
from rest_framework import viewsets
from rest_framework.decorators import detail_route
from rest_framework.response import Response

from .models import Edges, Meta, Nodes
from .serial import EdgesValueSerializer, MetaValueSerializer, TFValueSerializer


class MetaValueDistinctViewSet(viewsets.ReadOnlyModelViewSet):
    serializer_class = MetaValueSerializer
    queryset = Meta.objects.all()

    def list(self, request):
        queryset = Meta.objects.all().values("text").distinct()
        serializer = MetaValueSerializer(queryset, many=True)
        return Response(serializer.data)

    @detail_route()
    def searchType(self, request, pk=None):
        """"""
        uinput = request.query_params.get("uinput", None)
        queryset = []
        pk = pk.upper()
        if uinput != None:
            queryset = Meta.objects.filter(meta_type=pk,
                                           text__istartswith=uinput).all(

            ).values(
                "text").distinct()
        serializer = MetaValueSerializer(queryset, many=True)
        return Response(serializer.data)


class TFValueDistinctViewSet(viewsets.ReadOnlyModelViewSet):
    serializer_class = TFValueSerializer
    queryset = Nodes.objects.all().values("text").distinct()

    @detail_route()
    def searchName(self, request, pk=None):
        uinput = request.query_params.get("uinput", None)
        pk = pk.upper()
        if uinput and pk == "NODENAME":
            queryset = Nodes.objects.filter(
                Q(node_type="TF") | Q(node_type="-"),
                Q(text__istartswith=uinput)).all().values("text").distinct()
        elif uinput and pk == "NODETYPE":
            queryset = Nodes.objects.filter(
                Q(node_type="TF") | Q(node_type="-"),
                Q(text__istartswith=uinput)).all().values("text").distinct()
        else:
            queryset = Nodes.objects.filter(node_type="TF").all().values(
                "text").distinct()

        serializer = TFValueSerializer(queryset, many=True)

        return Response(serializer.data)


class EdgesValueDistinctViewSet(viewsets.ReadOnlyModelViewSet):
    serializer_class = EdgesValueSerializer
    queryset = Edges.objects.all().values("text").distinct()

    @detail_route()
    def searchName(self, request, pk=None):
        queryset = []
        uinput = request.query_params.get("uinput", None)
        pk = pk.upper()
        queryset = Edges.objects.filter(text__istartswith=uinput).all().values(
            "text").distinct()
        serializer = EdgesValueSerializer(queryset, many=True)
        return Response(serializer.data)
