from rest_framework import serializers

from querytgdb.models import TargetDBTF, Edges, MetaIddata
from .models import Meta, Nodes


class MetaValueSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = MetaIddata
        fields = ("meta_value",)


class TFValueSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = TargetDBTF
        fields = ("db_tf_agi",)


class EdgesValueSerializer(serializers.ModelSerializer):
    class Meta:
        model = Edges
        fields = ("edge_name",)


class TFTypeSerializer(serializers.HyperlinkedModelSerializer):
    text = serializers.CharField(source="node_type")

    class Meta:
        model = Nodes
        fields = ("text",)
