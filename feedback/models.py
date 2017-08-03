from typing import Type

from django.core.mail import send_mail
from django.db import models
from django.db.models.signals import post_save
from django.utils.translation import ugettext_lazy as _
from django.conf import settings


class Feedback(models.Model):
    """
    Records user feedback
    """
    name = models.CharField(max_length=255)
    feedback = models.TextField()
    created = models.DateTimeField(auto_now_add=True)


def send_feedback_alert(sender: Type[models.Model], instance: Feedback, *args, **kwargs) -> None:
    send_mail(
        _("%(name)s submitted feedback at %(time)s") % {
            'name': instance.name,
            'time': instance.created
        },
        instance.feedback,
        "noreply@coruzzilab-macpro.bio.nyu.edu",
        settings.ALERT_EMAILS
    )


post_save.connect(send_feedback_alert)
