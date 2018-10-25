from rest_framework import serializers

from querytgdb.models import Analysis as SavedAnalysis
from .models import Experiment


class AnalysisSerializer(serializers.Serializer):
    experiment = serializers.SlugField()
    analysis_id = serializers.SlugField()
    analysis_cutoff = serializers.CharField()
    analysis_method = serializers.CharField()
    analysis_batch = serializers.ChoiceField(('YES', 'NO'))
    analysis_notes = serializers.CharField()
    gene_list = serializers.FileField()
    experimental_design = serializers.FileField()


class ExperimentSerializer(serializers.Serializer):
    experiment_id = serializers.SlugField()
    analysis_id = serializers.SlugField()
    tf_id = serializers.SlugField()
    experiment = serializers.CharField()
    experiment_type = serializers.ChoiceField(("expression", "binding"))
    experiment_subtype = serializers.ChoiceField((
        "RNAseq",
        "Microarray",
        "4tU",
        "ChIPseq",
        "DamID"
    ))
    direction = serializers.ChoiceField((0, 1))
    genotype = serializers.CharField()
    data_source = serializers.CharField()
    time = serializers.IntegerField(min_value=0)
    growth_period = serializers.IntegerField(min_value=0)
    growth_medium = serializers.CharField()
    plasmid = serializers.CharField()
    control = serializers.CharField()
    tissue = serializers.ChoiceField(("shoot", "root"))
    treatments = serializers.CharField()
    replicates = serializers.CharField()
    batch = serializers.CharField()
    analysis_method = serializers.CharField()
    analysis_cutoff = serializers.CharField()
    analysis_command = serializers.CharField()
    analysis_batch = serializers.ChoiceField(("YES", "NO"))
    analysis_notes = serializers.CharField()
    tf_history_notes = serializers.CharField()
    experimenter = serializers.CharField()
    submission_date = serializers.DateField(input_formats=["%Y-%m-%d"])
    experiment_date = serializers.DateField(input_formats=["%Y-%m-%d"])
    metadata_notes = serializers.CharField()
    gene_list = serializers.FileField()
    expression_values = serializers.FileField()
    design = serializers.FileField()

    def validate_experiment_id(self, value):
        if Experiment.objects.filter(name=value).exists() or SavedAnalysis.objects.filter(name=value).exists():
            raise serializers.ValidationError("Experiment ID exstis!")

        return value
