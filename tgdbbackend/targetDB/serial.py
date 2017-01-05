from rest_framework import serializers

from .models import Edges, Meta, Nodes


class MetaValueSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = Meta
        fields = ["text", ]


class TFValueSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = Nodes
        fields = ["text", ]


class EdgesValueSerializer(serializers.ModelSerializer):
    class Meta:
        model = Edges
        fields = ["text"]


class TFTypeSerializer(serializers.HyperlinkedModelSerializer):
    text = serializers.CharField(source="node_type")

    class Meta:
        model = Nodes
        fields = ["text"]
