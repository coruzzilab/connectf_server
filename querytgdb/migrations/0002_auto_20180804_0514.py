# Generated by Django 2.0.7 on 2018-08-04 05:14

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('querytgdb', '0001_initial'),
    ]

    operations = [
        migrations.AlterField(
            model_name='analysisdata',
            name='value',
            field=models.CharField(max_length=200),
        ),
        migrations.AlterField(
            model_name='experimentdata',
            name='value',
            field=models.CharField(max_length=200),
        ),
    ]
