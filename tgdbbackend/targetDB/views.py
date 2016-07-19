from django.shortcuts import render

# Create your views here.
from rest_framework import viewsets
from .serial import MetaValueSerializer,TFValueSerializer,EdgesValueSerializer;
from .models import Meta,Nodes,Edges
from rest_framework.response import Response
from rest_framework.decorators import detail_route
from django.views.generic import View;

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
        queryset = [];
        pk = pk.upper();
        queryset =  Nodes.objects.filter(text__istartswith=pk).all().values("text").distinct();
        serializer = TFValueSerializer(queryset,many =True);
        return Response(serializer.data);
    
class EdgesValueDistinctViewSet(viewsets.ReadOnlyModelViewSet):
    serializer_class = EdgesValueSerializer;
    queryset = Edges.objects.all().values("text").distinct();
    
    @detail_route
    def searchName(self,request,pk=None):
        queryset = [];
        pk = pk.upper();
        queryset = Edges.objects.filter(text__istartswith=pk).all().values("text").distinct();
        serializer = EdgesValueSerializer(queryset,many=True);
        return Response(serializer.data);
    
class HandleQueryView(View):
    
    def get(self,request,*args,**kwargs):
        pass
    
    def post(self,request,*args,**kwargs):
        pass;