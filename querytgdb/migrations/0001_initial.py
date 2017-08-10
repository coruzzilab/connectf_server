# -*- coding: utf-8 -*-
# Generated by Django 1.11.3 on 2017-08-09 18:13
from __future__ import unicode_literals

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Analysis',
            fields=[
                ('analysis_p', models.AutoField(db_index=True, primary_key=True, serialize=False)),
                ('analysis_id', models.IntegerField(db_index=True, unique=True)),
                ('analysis_fullid', models.CharField(db_index=True, max_length=100)),
            ],
        ),
        migrations.CreateModel(
            name='AnalysisIddata',
            fields=[
                ('analysisid_id', models.AutoField(primary_key=True, serialize=False)),
                ('analysis_value', models.CharField(db_index=True, max_length=255)),
                ('analysis_type', models.CharField(db_index=True, max_length=100)),
                ('analysis_id', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='querytgdb.Analysis', to_field='analysis_id')),
            ],
        ),
        migrations.CreateModel(
            name='Annotation',
            fields=[
                ('ath_id', models.AutoField(db_index=True, primary_key=True, serialize=False)),
                ('agi_id', models.CharField(db_index=True, max_length=100)),
                ('ath_name', models.CharField(max_length=200)),
                ('ath_fullname', models.TextField(max_length=2000)),
                ('ath_gene_type', models.CharField(max_length=100)),
                ('ath_gene_fam', models.TextField(max_length=2000)),
            ],
        ),
        migrations.CreateModel(
            name='DAPdata',
            fields=[
                ('dap_interaction', models.AutoField(primary_key=True, serialize=False)),
                ('ath_id', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='querytgdb.Annotation')),
            ],
        ),
        migrations.CreateModel(
            name='Edges',
            fields=[
                ('edge_id', models.AutoField(db_index=True, primary_key=True, serialize=False)),
                ('edge_name', models.CharField(db_index=True, max_length=100, unique=True)),
            ],
        ),
        migrations.CreateModel(
            name='Interactions',
            fields=[
                ('interaction_id', models.AutoField(primary_key=True, serialize=False)),
            ],
        ),
        migrations.CreateModel(
            name='Metadata',
            fields=[
                ('meta_p', models.AutoField(db_index=True, primary_key=True, serialize=False)),
                ('meta_id', models.IntegerField(db_index=True, unique=True)),
                ('meta_fullid', models.CharField(db_index=True, max_length=100)),
            ],
        ),
        migrations.CreateModel(
            name='MetaIddata',
            fields=[
                ('metaid_id', models.AutoField(primary_key=True, serialize=False)),
                ('meta_value', models.CharField(db_index=True, max_length=255)),
                ('meta_type', models.CharField(db_index=True, max_length=100)),
                ('meta_id', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='querytgdb.Metadata', to_field='meta_id')),
            ],
        ),
        migrations.CreateModel(
            name='ReferenceId',
            fields=[
                ('ref_id', models.AutoField(primary_key=True, serialize=False)),
                ('analysis_id', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='querytgdb.Analysis', to_field='analysis_id')),
                ('meta_id', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='querytgdb.Metadata', to_field='meta_id')),
            ],
        ),
        migrations.CreateModel(
            name='Regulation',
            fields=[
                ('regulation_id', models.AutoField(primary_key=True, serialize=False)),
                ('foldchange', models.CharField(max_length=100)),
                ('pvalue', models.CharField(max_length=100)),
                ('ath_id', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='querytgdb.Annotation')),
                ('ref_id', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='querytgdb.ReferenceId')),
            ],
        ),
        migrations.CreateModel(
            name='TargetDBTF',
            fields=[
                ('db_tf_id', models.AutoField(db_index=True, primary_key=True, serialize=False)),
                ('db_tf_agi', models.CharField(db_index=True, max_length=100, unique=True)),
            ],
        ),
        migrations.AddField(
            model_name='interactions',
            name='db_tf_id',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='querytgdb.TargetDBTF'),
        ),
        migrations.AddField(
            model_name='interactions',
            name='edge_id',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='querytgdb.Edges'),
        ),
        migrations.AddField(
            model_name='interactions',
            name='ref_id',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='querytgdb.ReferenceId'),
        ),
        migrations.AddField(
            model_name='interactions',
            name='target_id',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='querytgdb.Annotation'),
        ),
        migrations.AddField(
            model_name='dapdata',
            name='db_tfid',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='querytgdb.TargetDBTF'),
        ),
    ]
