from rest_framework import serializers

from querytgdb.models import Edges, MetaIddata, TargetDBTF
# from .models import Nodes


class MetaValueSerializer(serializers.ModelSerializer):
    value = serializers.CharField(source='meta_value')

    class Meta:
        model = MetaIddata
        fields = ("value",)


class TFValueSerializer(serializers.ModelSerializer):
    value = serializers.CharField(source='db_tf_agi')
    name = serializers.CharField(required=False, source='ath_name')

    class Meta:
        model = TargetDBTF
        fields = ("value", "name")


class EdgesValueSerializer(serializers.ModelSerializer):
    value = serializers.CharField(source='edge_name')

    class Meta:
        model = Edges
        fields = ("value",)


# class TFTypeSerializer(serializers.ModelSerializer):
#     text = serializers.CharField(source="node_type")
#
#     class Meta:
#         model = Nodes
#         fields = ("text",)
