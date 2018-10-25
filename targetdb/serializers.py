from rest_framework import serializers

from querytgdb.models import Analysis


class TFValueSerializer(serializers.ModelSerializer):
    value = serializers.CharField(source='gene_id')
    name = serializers.CharField(required=False, source='gene_name')

    class Meta:
        model = Analysis
        fields = ("value", "name")
