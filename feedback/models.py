from django.db import models


class Feedback(models.Model):
    """
    Records user feedback
    """
    name = models.CharField(max_length=255)
    feedback = models.TextField()
    created = models.DateTimeField(auto_now_add=True)
