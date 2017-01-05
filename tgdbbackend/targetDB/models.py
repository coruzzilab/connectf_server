from django.db import models


# Create your models here.
class Meta(models.Model):
    meta_p = models.IntegerField()
    meta_id = models.CharField(max_length=100, primary_key=True)
    text = models.CharField(max_length=100, db_column="meta_value")
    meta_type = models.CharField(max_length=100)

    class Meta:
        verbose_name = "Meta"
        verbose_name_plural = "Metas"
        db_table = "Meta"
        managed = False


class Nodes(models.Model):
    node_id = models.IntegerField()
    text = models.CharField(max_length=100, db_column="node_name")
    node_type = models.CharField(max_length=100, db_column="node_type")

    class Meta:
        verbose_name = "Node"
        verbose_name_plural = "Nodes"
        db_table = "Nodes"
        managed = False


class Edges(models.Model):
    edge_id = models.IntegerField()
    text = models.CharField(max_length=100, db_column="edge_name")

    class Meta:
        verbose_name = "Edge"
        verbose_name_plural = "Edges"
        db_table = "Edges"
        managed = False


class Interactions(models.Model):
    interaction_id = models.IntegerField()
    node_1_id = models.IntegerField()
    node_2_id = models.IntegerField()
    edge_id = models.IntegerField()
    meta_id = models.CharField(max_length=100)

    class Meta:
        verbose_name = "Interaction"
        verbose_name_plural = "Interactions"
        managed = False
        db_table = "Interactions"
