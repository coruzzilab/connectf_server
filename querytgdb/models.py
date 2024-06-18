import re
import sys
from collections import defaultdict
from typing import DefaultDict

import numpy as np
from django.db import models
from django.utils.text import gettext_lazy as _


class Analysis(models.Model):
    """
    Individual analyses
    """
    tf = models.ForeignKey("Annotation", on_delete=models.CASCADE)

    @property
    def meta_dict(self) -> DefaultDict[str, str]:
        return defaultdict(str, self.analysisdata_set.values_list('key__name', 'value'))

    @property
    def name(self) -> str:
        return '{0.tech}_{0.analysis_method}_{0.pk}'.format(self)

    @property
    def tech(self) -> str:
        d = self.meta_dict

        return d['EXPERIMENTER'] + d['DATE'].replace('-', '') + '_' + d['TECHNOLOGY']

    @property
    def analysis_method(self) -> str:
        d = self.meta_dict

        return d['ANALYSIS_METHOD'] + '_' + re.sub(r'\s+', '', d['ANALYSIS_CUTOFF'])

    class Meta:
        verbose_name_plural = "analyses"


class MetaKey(models.Model):
    """
    Keys for metadata fields
    """
    name = models.CharField(max_length=100, db_index=True, unique=True)
    searchable = models.BooleanField(default=True)

    def __str__(self):
        return "{0.name} searchable={0.searchable}".format(self)


class AnalysisData(models.Model):
    """
    Metadata for an analysis
    """
    analysis = models.ForeignKey(Analysis, on_delete=models.CASCADE)
    key = models.ForeignKey(MetaKey, on_delete=models.CASCADE)
    value = models.CharField(max_length=200)

    class Meta:
        unique_together = (('analysis', 'key'),)


class Annotation(models.Model):
    """
    Gene annotations
    """
    gene_id = models.CharField(max_length=100, db_index=True, unique=True)
    name = models.CharField(max_length=100)
    fullname = models.TextField(max_length=2000)
    gene_type = models.CharField(max_length=100)
    gene_family = models.TextField(max_length=2000)

    def __str__(self):
        return self.name or self.gene_id

    @property
    def gene_name_symbol(self) -> str:
        return self.gene_id + (f" ({self.name})" if self.name else "")


class EdgeData(models.Model):
    """
    Extra edge property relationships (Unrelated to experiments)
    """
    tf = models.ForeignKey(Annotation, on_delete=models.CASCADE, related_name='tfs')
    target = models.ForeignKey(Annotation, on_delete=models.CASCADE, related_name='targets')
    type = models.ForeignKey("EdgeType", on_delete=models.CASCADE)

    class Meta:
        unique_together = (("tf", "target", "type"),)


class EdgeType(models.Model):
    """
    Extra edge property types (DAP, DAPamp, etc.)
    """
    name = models.CharField(max_length=100, unique=True)
    directional = models.BooleanField(default=True)


class Interaction(models.Model):
    """
    Edge for TF—target interactions
    """
    analysis = models.ForeignKey(Analysis, on_delete=models.CASCADE)
    target = models.ForeignKey(Annotation, on_delete=models.CASCADE)

    class Meta:
        unique_together = (('analysis', 'target'),)


class NumpyFloat(models.FloatField):
    description = _("Float field that handles numpy types.")

    def from_db_value(self, value, expression, connection):
        if value is None:
            return np.nan
        return value

    def to_python(self, value):
        if np.isnan(value):
            return np.nan

        return super().to_python(value)

    def get_prep_value(self, value):
        value = super().get_prep_value(value)

        if value is None or np.isnan(value):
            return None
        elif abs(value) == sys.float_info.max:
            # The mysql connector makes some odd rounding errors, so convert to string before passing it along the way.
            return str(value)

        return value


class Regulation(models.Model):
    """
    P-value and fold change for TF—target interactions

    Not necessary for every TF analysis
    """
    analysis = models.ForeignKey(Analysis, on_delete=models.CASCADE)
    target = models.ForeignKey(Annotation, on_delete=models.CASCADE)
    foldchange = NumpyFloat(null=True)
    p_value = models.FloatField()

    class Meta:
        unique_together = (('analysis', 'target'),)

class ImportHistory(models.Model):
    """
    History of imported data
    """
    fileName = models.FilePathField()
    type = models.CharField(max_length=10)
    imported_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return f"{self.type}: {self.fileName} imported at {self.imported_at}"