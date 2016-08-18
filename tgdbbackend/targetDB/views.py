from django.shortcuts import render

# Create your views here.
from rest_framework import viewsets
from .serial import MetaValueSerializer,TFValueSerializer,EdgesValueSerializer,TFTypeSerializer;
from .models import Meta,Nodes,Edges
from rest_framework.response import Response
from rest_framework.decorators import detail_route
from django.views.generic import View;
from django.db.models import Q

class MetaValueDistinctViewSet(viewsets.ReadOnlyModelViewSet):
    serializer_class = MetaValueSerializer;
    queryset = Meta.objects.all();
    
    def list(self, request):
        queryset = Meta.objects.all().values("text").distinct();
        serializer = MetaValueSerializer(queryset, many=True);
        return Response(serializer.data);
    
    @detail_route()
    def searchType(self, request,pk=None):
        """"""
        uinput = request.query_params.get("uinput",None);
        queryset=[];
        pk = pk.upper();
        if uinput!=None:
            queryset = Meta.objects.filter(meta_type=pk,text__istartswith=uinput).all().values("text").distinct();
        serializer=MetaValueSerializer(queryset,many=True);
        return Response(serializer.data);
    
class TFValueDistinctViewSet(viewsets.ReadOnlyModelViewSet):
    serializer_class = TFValueSerializer;
    queryset = Nodes.objects.all().values("text").distinct();
    
    @detail_route()
    def searchName(self,request,pk=None):
        uinput = request.query_params.get("uinput",None);
        queryset = [];
        pk = pk.upper();
        if uinput!=None and pk=="NODENAME":
            queryset =  Nodes.objects.filter(Q(node_type="TF") | Q(node_type="-"),Q(text__istartswith=uinput)).all().values("text").distinct();
            serializer = TFValueSerializer(queryset,many =True);
        elif uinput!=None and pk=="NODETYPE":
            queryset = Nodes.objects.filter(Q(node_type=="TF") | Q(node_type="-"),Q(text__istartswith=uinput)).all().values("text").distinct();   
            serializer = TFValueSerializer(queryset,many =True);
                 
        return Response(serializer.data);
    
class EdgesValueDistinctViewSet(viewsets.ReadOnlyModelViewSet):
    serializer_class = EdgesValueSerializer;
    queryset = Edges.objects.all().values("text").distinct();
    
    @detail_route()
    def searchName(self,request,pk=None):
        queryset = [];
        uinput = request.query_params.get("uinput",None);
        pk = pk.upper();
        queryset = Edges.objects.filter(text__istartswith=uinput).all().values("text").distinct();
        serializer = EdgesValueSerializer(queryset,many =True);
        return Response(serializer.data);
    
