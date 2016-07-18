from rest_framework import serializers;

from .models import Meta,Nodes,Edges,Interactions;

class MetaValueSerializer(serializers.HyperlinkedModelSerializer):
    
    class Meta:
        model = Meta;
        fields = ["text",]
        