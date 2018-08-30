from django.db import models


class Edge(models.Model):
    """
    Target Edge names
    """
    name = models.CharField(max_length=100, unique=True)

    def __str__(self):
        return "{} ({})".format(self.name, self.id)


class Experiment(models.Model):
    """
    A set of analyses for a TF
    """
    name = models.CharField(max_length=100, unique=True)
    tf = models.ForeignKey("Annotation", on_delete=models.CASCADE)

    def __str__(self):
        return "{} ({})".format(self.name, self.id)


class Analysis(models.Model):
    """
    Individual analyses
    """
    name = models.CharField(max_length=100, db_index=True)
    experiment = models.ForeignKey(Experiment, on_delete=models.CASCADE)

    def __str__(self):
        return "{} ({})".format(self.name, self.id)

    class Meta:
        unique_together = (('experiment', 'name'),)


class ExperimentData(models.Model):
    """
    Metadata for an Experiment
    """
    experiment = models.ForeignKey(Experiment, on_delete=models.CASCADE)
    key = models.CharField(max_length=100, db_index=True)
    value = models.CharField(max_length=200)

    class Meta:
        unique_together = (('experiment', 'key'),)


class AnalysisData(models.Model):
    """
    Metadata for an analysis
    """
    analysis = models.ForeignKey(Analysis, on_delete=models.CASCADE)
    key = models.CharField(max_length=100, db_index=True)
    value = models.CharField(max_length=200)

    class Meta:
        unique_together = (('analysis', 'key'),)


class Annotation(models.Model):
    """
    Gene annotations
    """
    gene_id = models.CharField(max_length=100, unique=True)
    name = models.CharField(max_length=100)
    fullname = models.TextField(max_length=2000)
    gene_type = models.CharField(max_length=100)
    gene_family = models.TextField(max_length=2000)

    def __str__(self):
        return self.name or self.gene_id


class EdgeData(models.Model):
    """
    Extra edge property relationships
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


class Interaction(models.Model):
    """
    Edge for TF—target interactions
    """
    analysis = models.ForeignKey(Analysis, on_delete=models.CASCADE)
    target = models.ForeignKey(Annotation, on_delete=models.CASCADE)
    edge = models.ForeignKey(Edge, on_delete=models.CASCADE)

    class Meta:
        unique_together = (('analysis', 'target'),)


class Regulation(models.Model):
    """
    P-value and fold change for TF—target interactions

    Not necessary for every TF analysis
    """
    analysis = models.ForeignKey(Analysis, on_delete=models.CASCADE)
    target = models.ForeignKey(Annotation, on_delete=models.CASCADE)
    foldchange = models.FloatField()
    p_value = models.FloatField()

    class Meta:
        unique_together = (('analysis', 'target'),)
