# Generated by Django 2.0.7 on 2018-07-24 15:47

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('upload', '0001_initial'),
    ]

    operations = [
        migrations.AlterField(
            model_name='experiment',
            name='name',
            field=models.CharField(db_index=True, max_length=128, unique=True),
        ),
    ]
