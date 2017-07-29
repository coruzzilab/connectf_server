import io
import uuid
from datetime import date

from django import forms
from django.core.files.storage import FileSystemStorage
from django.core.mail import send_mail
from django.utils.translation import ugettext_lazy as _

storage = FileSystemStorage()

META_DATA = """Experiment_ID: {experiment_id}
Transcription_Factor_ID: {tf_id}
Experiment: {experiment}
Experiment_Type: {experiment_type}
Expression_Type: {expression_type}
Binding_Type: {binding_type}
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


class UploadForm(forms.Form):
    experiment_id = forms.SlugField(widget=forms.TextInput(attrs={
        'placeholder': 'e.g. AT4G24020_AS090116_RNASEQ ('
                       'TFID_ExperimenterInitials&ExperimentDate_Type)'
    }))

    tf_id = forms.CharField(widget=forms.TextInput(attrs={
        'placeholder': 'e.g. AT4G24020'
    }))

    experiment = forms.CharField(widget=forms.TextInput(attrs={
        'placeholder': 'e.g. Target/Inplanta etc.'
    }))

    experiment_type = forms.ChoiceField(choices=[
        ("Expression", "Expression"),
        ("Binding", "Binding")
    ])

    expression_type = forms.ChoiceField(choices=[
        ("NA", "NA"),
        ("RNAseq", "RNAseq"),
        ("Microarray", "Microarray"),
        ("4tU", "4tU")
    ])

    binding_type = forms.ChoiceField(choices=[
        ("NA", "NA"),
        ("ChIPseq", "ChIPseq"),
        ("DamID", "DamID")
    ])

    direction = forms.CharField(widget=forms.TextInput(attrs={
        'placeholder': 'e.g. 1'
    }))

    genotype = forms.CharField(widget=forms.TextInput(attrs={
        'placeholder': 'e.g. nlp7-1 or Col-0'
    }))

    data_source = forms.CharField(widget=forms.TextInput(attrs={
        'placeholder': 'e.g. CoruzziLab_unpublished_AS ('
                       'AS=ExperimenterInitials)'
    }))

    time = forms.IntegerField(
        min_value=0,
        widget=forms.NumberInput(attrs={'placeholder': 0})
    )

    growth_period = forms.IntegerField(
        min_value=0,
        widget=forms.NumberInput(attrs={'placeholder': 0})
    )

    growth_medium = forms.CharField(widget=forms.TextInput(attrs={
        'placeholder': 'e.g. plates,1mM_KNO3'
    }))

    plasmid = forms.CharField(widget=forms.TextInput(attrs={
        'placeholder': 'e.g. NLP7pBOB11'
    }))

    control = forms.ChoiceField(choices=[
        ("mDEX", "mDEX"),
        ("EmptyVector", "EmptyVector")
    ])

    treatments = forms.CharField(widget=forms.TextInput(attrs={
        'placeholder': 'e.g. PN,MDEX,PDEX,PCHX'
    }))

    replicates = forms.CharField(widget=forms.TextInput(attrs={
        'placeholder': 'e.g. PN+PCHX(3),PN+PDEX+PCHX(3)'
    }))

    batch = forms.CharField(widget=forms.TextInput(attrs={
        'placeholder': 'List the exp ids for experiments done together (comma '
                       'separated)',
        'style': 'width: 400px;'
    }))

    analysis_method = forms.CharField(widget=forms.TextInput(attrs={
        'placeholder': 'e.g. DESEQ2'
    }))
    analysis_cutoff = forms.CharField(widget=forms.TextInput(attrs={
        'placeholder': 'e.g. FDR<0.1'
    }))
    analysis_command = forms.CharField(widget=forms.TextInput(attrs={
        'placeholder': 'e.g. aov('
                       'dataframe$y~dataframe$Nitrogen*dataframe$Genotype'
                       '*dataframe$Tissue)'
    }))
    analysis_notes = forms.CharField(widget=forms.Textarea(attrs={
        'placeholder': 'e.g. Add notes about the analysis (alignment method, '
                       'annotation version, read count tool  etc.)'
    }))
    tf_history_notes = forms.CharField(widget=forms.Textarea(attrs={
        'placeholder': 'e.g. Wang_2009_Plant_physiology, '
                       'Castaings_2009_Plant_journal'
    }))

    experimenter = forms.ChoiceField(choices=[
        ("Anna_Schinke", "Anna Schinke"),
        ("Matthew_Brooks", "Matthew Brooks"),
        ("Ying_Li", "Ying Li"),
        ("Joseph_Swift", "Joseph Swift"),
        ("Sophie_Leran", "Sophie Leran"),
        ("Joan_Doidy", "Joan Doidy"),
        ("Eleonore_Bouguyon", "Eleonore Bouguyon"),
        ("Chia-Yi_Cheng", "Chia-Yi Cheng")
    ])

    submission_date = forms.DateField(
        initial=date.today,
        input_formats=["%m-%d-%Y"],
        widget=forms.DateInput(format="%m-%d-%Y",
                               attrs={'placeholder': 'mm-dd-yyyy'})
    )
    experiment_date = forms.DateField(
        input_formats=["%m-%d-%Y"],
        widget=forms.DateInput(attrs={'placeholder': 'mm-dd-yyyy'})
    )

    metadata_notes = forms.CharField(widget=forms.Textarea(attrs={
        'placeholder': 'e.g. Add notes for this data'
    }))

    gene_list = forms.FileField()

    expression_values = forms.FileField()

    design = forms.FileField()

    def save_files(self):
        uid = uuid.uuid1()
        self.save_metadata(uid)
        self.handle_file(self.cleaned_data['gene_list'], "genelist", uid)
        self.handle_file(self.cleaned_data['expression_values'], "rawdata", uid)
        self.handle_file(self.cleaned_data['design'], "expdesign", uid)

    def handle_file(self, upload_file, suffix, uid):
        storage.save(
            '{experiment_id}_{uuid}_{suffix}.txt'.format(
                uuid=uid,
                suffix=suffix,
                experiment_id=self.cleaned_data['experiment_id']),
            upload_file)

    def save_metadata(self, uid):
        storage.save('{experiment_id}_{uid}_metadata.txt'.format(
            experiment_id=self.cleaned_data['experiment_id'],
            uid=uid),
            io.StringIO(META_DATA.format(**self.cleaned_data))
        )

    def send_mail(self):
        experimenter = self.cleaned_data['experimenter']
        send_mail(
            _("%(name)s has uploaded an experiment") % {"name": experimenter},
            "as titled",
            "noreply@coruzzilab-macpro.bio.nyu.edu",
            [
                "rt76@nyu.edu",
                "clj327@nyu.edu"
            ]
        )
