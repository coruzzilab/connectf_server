from rest_framework import serializers

from querytgdb.models import Annotation, Edges, MetaIddata, Metadata, TargetDBTF


class MetaValueSerializer(serializers.ModelSerializer):
    value = serializers.CharField(source='meta_value')

    class Meta:
        model = MetaIddata
        fields = ("value",)


class ExperimentIdSerializer(serializers.ModelSerializer):
    value = serializers.CharField(source='meta_fullid')

    class Meta:
        model = Metadata
        fields = ("value",)


class AnnotationSerializer(serializers.ModelSerializer):
    value = serializers.CharField(source='agi_id')
    name = serializers.CharField(required=False, source='ath_name')

    class Meta:
        model = Annotation
        fields = ("value", "name")


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
