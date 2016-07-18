from django.shortcuts import render

# Create your views here.
from rest_framework import viewsets
from .serial import MetaValueSerializer;
from .models import Meta
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
        if uinput!=None:
            queryset = Meta.objects.filter(meta_type=pk,text__istartswith=uinput).all().values("text").distinct();
        serializer=MetaValueSerializer(queryset,many=True);
        return Response(serializer.data);
    

class HandleQueryView(View):
    
    def get(self,request,*args,**kwargs):
        pass;   
    
    def post(self,request,*args,**kwargs):
        pass;