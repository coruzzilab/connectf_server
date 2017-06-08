#!/usr/bin/python
# -*- coding: utf-8 -*-

'''

This script setup the Django mysql database for TargetDB version2

'''

from django.db import models


class TargetDBTF(models.Model):

    db_tf_id = models.AutoField(primary_key=True, db_index=True)
    db_tf_agi = models.CharField(max_length=100, db_index=True, unique=True)


class Edges(models.Model):

    edge_id = models.AutoField(primary_key=True, db_index=True)
    edge_name = models.CharField(max_length=100, db_index=True, unique=True)


class Metadata(models.Model):

    meta_p= models.AutoField(primary_key=True, db_index=True)
    # one-to-many relationship: one metaid can have multiple referenceids
    meta_id = models.IntegerField(db_index = True, unique= True)
    meta_fullid= models.CharField(max_length=100, db_index=True)


class Analysis(models.Model):

    analysis_p= models.AutoField(primary_key=True, db_index=True)
    # one-to-one relationship: one analysisid will have one referenceid
    analysis_id= models.IntegerField(db_index=True, unique= True)
    analysis_fullid= models.CharField(max_length=100, db_index=True)


class ReferenceId(models.Model):

    ref_id= models.AutoField(primary_key=True)
    meta_id = models.ForeignKey(Metadata, on_delete=models.CASCADE, to_field= 'meta_id')
    analysis_id = models.ForeignKey(Analysis, on_delete=models.CASCADE, to_field= 'analysis_id')


class MetaIddata(models.Model):

    metaid_id= models.AutoField(primary_key=True)
    meta_value= models.CharField(max_length=255, db_index=True)
    meta_type= models.CharField(max_length=100, db_index=True)
    meta_id = models.ForeignKey(Metadata, on_delete=models.CASCADE, to_field='meta_id')


class AnalysisIddata(models.Model):

    analysisid_id= models.AutoField(primary_key=True)
    analysis_value=  models.CharField(max_length=255, db_index=True)
    analysis_type= models.CharField(max_length=100, db_index=True)
    analysis_id= models.ForeignKey(Analysis, on_delete=models.CASCADE, to_field= 'analysis_id')


class Annotation(models.Model):

    ath_id= models.AutoField(primary_key=True, db_index=True)
    agi_id= models.CharField(max_length=100, db_index=True)
    ath_name= models.CharField(max_length=200)
    ath_fullname= models.TextField(max_length=2000)
    ath_gene_type= models.CharField(max_length=100)
    ath_gene_fam= models.TextField(max_length=2000)


class DAPdata(models.Model):

    dap_interaction= models.AutoField(primary_key= True)
    # one-to-many relationship: one TF can have multiple interactions
    db_tfid= models.ForeignKey(TargetDBTF, on_delete=models.CASCADE, to_field= 'db_tf_id')
    # one-to-many relationship: one target can have multiple interactionids
    ath_id= models.ForeignKey(Annotation, on_delete=models.CASCADE, to_field= 'ath_id')


class Interactions(models.Model):

    interaction_id= models.AutoField(primary_key=True)
    # one-to-many relationship: one TF can have multiple interactions
    db_tf_id= models.ForeignKey(TargetDBTF, on_delete=models.CASCADE, to_field= 'db_tf_id')
    # one-to-many relationship: one target can have multiple interactionids
    target_id= models.ForeignKey(Annotation, on_delete=models.CASCADE, to_field= 'ath_id')
    # one-to-many relationship: one edge can have multiple interactionid
    edge_id= models.ForeignKey(Edges, on_delete=models.CASCADE, to_field= 'edge_id')
    # one-to-many relationship: one referenceid can have multiple interactionid
    ref_id=models.ForeignKey(ReferenceId, on_delete=models.CASCADE, to_field= 'ref_id')


class Regulation(models.Model):

    regulation_id= models.AutoField(primary_key=True)
    # one-to-many relationship: one referenceid can have multiple regulationids
    ref_id= models.ForeignKey(ReferenceId, on_delete=models.CASCADE, to_field= 'ref_id')
    # one-to-many relationship:one AGI ID can have multiple regulationids
    ath_id= models.ForeignKey(Annotation, on_delete=models.CASCADE, to_field= 'ath_id')
    foldchange= models.CharField(max_length=100)
    pvalue= models.CharField(max_length=100)

