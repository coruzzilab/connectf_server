from __future__ import unicode_literals

from django.db import models

# Create your models here.
class Meta(models.Model):
    meta_p =  models.IntegerField();
    meta_id =  models.CharField(max_length=100,primary_key=True);
    text = models.CharField(max_length=100,db_column = "meta_value");
    meta_type = models.CharField(max_length=100);
    
    class Meta:
        verbose_name = "Meta"
        verbose_name_plural = "Metas"
        db_table = "Meta";
        managed = False;
    def __str__(self):
        pass

class Nodes(models.Model):
    node_id = models.IntegerField();
    text = models.CharField(max_length=100, db_column = "node_name");
    node_type = models.CharField(max_length=100);

    class Meta:
        verbose_name = "Nodes"
        verbose_name_plural = "Nodess"
	db_table = "Nodes"
	managed = False;
    def __str__(self):
        pass
    
class Edges(models.Model):
    edge_id = models.IntegerField();
    text = models.CharField(max_length=100, db_column = "edge_name");

    class Meta:
        verbose_name = "Edges"
        verbose_name_plural = "Edgess"
	db_table = "Edges"
	managed = False;
    def __str__(self):
        pass
    
class Interactions(models.Model):
    interaction_id = models.IntegerField();
    node_1_id = models.IntegerField();
    node_2_id = models.IntegerField();
    edge_id = models.IntegerField();
    meta_id = models.CharField(max_length=100);
    class Meta:
        verbose_name = "Interactions"
        verbose_name_plural = "Interactionss"
	managed = False;
    def __str__(self):
        pass
    