from django.db import models


# Create your models here.
class Experiment(models.Model):
    name = models.CharField(max_length=128, unique=True, db_index=True)

    def __str__(self):
        return self.name
