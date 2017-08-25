import io
import uuid

from django.core.files.storage import FileSystemStorage
from django.views.generic import FormView
from rest_framework.generics import CreateAPIView
from rest_framework.permissions import AllowAny
from rest_framework.serializers import ValidationError

from querytgdb.utils import validate_autosubmit as validate
from querytgdb.utils.insert_tgdb import insertdata
from . import forms, models, serializers

storage = FileSystemStorage()

META_DATA = """\
Experiment_ID: {experiment_id}
Transcription_Factor_ID: {tf_id}
Experiment: {experiment}
Experiment_Type: {experiment_type}
Experiment_Subtype: {experiment_subtype}
Direction: {direction}
Genotype: {genotype}
Data_Source: {data_source}
Time: {time}
Growth_Period: {growth_period}
Growth_Medium: {growth_medium}
Plasmids: {plasmid}
Control: {control}
Treatments: {treatments}
Replicates: {replicates}
Batch: {batch}
Analysis_method: {analysis_method}
Analysis_cutoff: {analysis_cutoff}
Analysis_command: {analysis_command}
Analysis_batch: {analysis_batch}
Analysis_notes: {analysis_notes}
TF_History: {tf_history_notes}
Experimenter: {experimenter}
Submission_date: {submission_date:%m-%d-%Y}
Experiment_date: {experiment_date:%m-%d-%Y}
Metadata_Notes: {metadata_notes}
file: {gene_list.name}
file: {expression_values.name}
file: {design.name}
"""

META_DATA_ANALYSIS = """\
Experiment_ID: {experiment}
Analysis_ID: {analysis_id}
Analysis_method: {analysis_method}
Analysis_cutoff: {analysis_cutoff}
Analysis_batch: {analysis_batch}
Analysis_notes: {analysis_notes}
file: {gene_list.name}
file: {experimental_design.name}
"""


def save_metadata(uid, data, meta_str=META_DATA):
    name = storage.save('{experiment_id}_{uid}_metadata.txt'.format(
        experiment_id=data['experiment_id'],
        uid=uid),
        io.StringIO(meta_str.format(**data))
    )

    return storage.path(name)


def handle_file(upload_file, data, suffix, uid):
    name = storage.save(
        '{experiment_id}_{uuid}_{suffix}.txt'.format(
            uuid=uid,
            suffix=suffix,
            experiment_id=data['experiment_id']),
        upload_file)

    return storage.path(name)


def handle_analysis_file(validated_data, key, suffix, uid):
    storage.save(
        '{experiment}_{analysis_id}_{uuid}_{suffix}.txt'.format(
            uuid=uid,
            suffix=suffix,
            experiment=validated_data['experiment'],
            analysis_id=validated_data['analysis_id']),
        validated_data[key])


def save_analysis_metadata(validated_data, uid, meta_str=META_DATA_ANALYSIS):
    storage.save('{experiment}_{analysis_id}_{uid}_metadata.txt'.format(
        experiment=validated_data['experiment'],
        analysis_id=validated_data['analysis_id'],
        uid=uid),
        io.StringIO(meta_str.format(**validated_data))
    )


# Create your views here.
class UploadView(FormView):
    template_name = 'upload/index.html'
    form_class = forms.ExperimentUploadForm
    success_url = 'http://coruzzilab-macpro.bio.nyu.edu'

    def form_valid(self, form: forms.ExperimentUploadForm):
        form.save_files()
        form.send_mail()
        return super().form_valid(form)


class UploadAnalysisView(CreateAPIView):
    serializer_class = serializers.AnalysisSerializer
    permission_classes = (AllowAny,)

    def perform_create(self, serializer):
        data = serializer.validated_data

        uid = uuid.uuid1()

        save_analysis_metadata(data, uid)
        handle_analysis_file(data, 'gene_list', "genelist", uid)
        handle_analysis_file(data, 'experimental_design', "expdesign", uid)


class UploadExperimentView(CreateAPIView):
    serializer_class = serializers.ExperimentSerializer
    permission_classes = (AllowAny,)

    def perform_create(self, serializer):
        data = serializer.validated_data

        uid = uuid.uuid1()
        try:
            # save data even if checks don't pass
            meta_path = save_metadata(uid, data)
            gene_list_path = handle_file(data['gene_list'], data, "genelist", uid)
            exp_value_path = handle_file(data['expression_values'], data, "rawdata", uid)
            design_path = handle_file(data['design'], data, "expdesign", uid)
            models.Experiment.objects.create(name=data['experiment_id'])

            # actually check
            metadict = validate.validate_metadata(meta_path)
            validate.validate_genelist(gene_list_path, metadict)
            validate.validate_readcount_expdesign(exp_value_path, design_path)

            # insert if checks pass
            insertdata(gene_list_path, meta_path,
                       '/Users/Reetu/Documents/Projects/TargetDB_V2/170801_TargetDB_latestdata/TargetDBdata/dap-seq'
                       '.all.txt')

        except ValueError as e:
            raise ValidationError(str(e)) from e
