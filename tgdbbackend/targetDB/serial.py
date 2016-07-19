from rest_framework import serializers;

from .models import Meta,Nodes,Edges,Interactions;

class MetaValueSerializer(serializers.HyperlinkedModelSerializer):
    
    class Meta:
        model = Meta;
        fields = ["text",]

class TFValueSerializer(serializers.HyperlinkedModelSerializer):
    
    class Meta:
        model = Nodes
        fields = ["text",]
        
class EdgesValueSerializer(serializers.HyperlinkedModelSerializer):
    
    class Meta:
        model = Edges
        field = ["text",]