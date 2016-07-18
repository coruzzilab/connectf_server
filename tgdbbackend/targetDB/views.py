from django.shortcuts import render

# Create your views here.
from rest_framework import viewsets
from .serial import MetaValueSerializer;
from .models import Meta
from rest_framework.response import Response
from rest_framework.decorators import detail_route

class MetaValueDistinctViewSet(viewsets.ReadOnlyModelViewSet):
    serializer_class = MetaValueSerializer;
    queryset = Meta.objects.all();
    def list(self, request):
        queryset = Meta.objects.all().values("text").distinct();
        serializer = MetaValueSerializer(queryset, many=True);
        return Response(serializer.data);
    
    def searchType(self, request,area=None, uinput=""):
        """"""
        queryset = Meta.objects.filter(meta_type=area,text__startswith=uinput).all().values("text").distinct();
        serializer=MetaValueSerializer(queryset);
        return Response(serializer.data);
    
