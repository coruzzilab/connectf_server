import io
import uuid

from rest_framework import serializers

from .forms import storage

META_DATA = """\
Experiment_ID: {experiment}
Analysis_ID: {analysis_id}
Analysis_method: {analysis_method}
Analysis_cutoff: {analysis_cutoff}
Analysis_batch: {analysis_batch}
Analysis_notes: {analysis_notes}
file: {gene_list.name}
file: {experimental_design.name}
"""


class AnalysisSerializer(serializers.Serializer):
    experiment = serializers.SlugField()
    analysis_id = serializers.SlugField()
    analysis_cutoff = serializers.CharField()
    analysis_method = serializers.CharField()
    analysis_batch = serializers.CharField()
    analysis_notes = serializers.CharField()
    gene_list = serializers.FileField()
    experimental_design = serializers.FileField()

    def create(self, validated_data):
        self.save_files(validated_data)
        return validated_data

    def save_files(self, validated_data):
        uid = uuid.uuid1()
        self.save_metadata(validated_data, uid)
        self.handle_file(validated_data, 'gene_list', "genelist", uid)
        self.handle_file(validated_data, 'experimental_design', "expdesign", uid)

    def handle_file(self, validated_data, key, suffix, uid):
        storage.save(
            '{experiment}_{analysis_id}_{uuid}_{suffix}.txt'.format(
                uuid=uid,
                suffix=suffix,
                experiment=validated_data['experiment'],
                analysis_id=validated_data['analysis_id']),
            validated_data[key])

    def save_metadata(self, validated_data, uid, meta_str=META_DATA):
        storage.save('{experiment}_{analysis_id}_{uid}_metadata.txt'.format(
            experiment=validated_data['experiment'],
            analysis_id=validated_data['analysis_id'],
            uid=uid),
            io.StringIO(meta_str.format(**validated_data))
        )
