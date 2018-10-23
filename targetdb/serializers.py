from rest_framework import serializers

from querytgdb.models import Experiment


class TFValueSerializer(serializers.ModelSerializer):
    value = serializers.CharField(source='gene_id')
    name = serializers.CharField(required=False, source='gene_name')

    class Meta:
        model = Experiment
        fields = ("value", "name")
