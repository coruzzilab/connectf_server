# Generated by Django 2.1.5 on 2019-01-17 18:50

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('querytgdb', '0006_remove_analysis_name'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='analysis',
            options={'verbose_name_plural': 'analyses'},
        ),
        migrations.AddField(
            model_name='edgetype',
            name='directional',
            field=models.BooleanField(default=True),
        ),
    ]
