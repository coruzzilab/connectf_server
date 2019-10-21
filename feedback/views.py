import logging

import boto3
import requests
from botocore.exceptions import ClientError
from django.conf import settings
from django.forms import CharField
from django.forms.models import ModelForm, ValidationError
from django.http import HttpRequest, HttpResponse, HttpResponseBadRequest
from django.views import View

from . import models

logger = logging.getLogger(__name__)

CHARSET = 'UTF-8'


def aws_send_email(instance: models.Feedback):
    client = boto3.client('sns',
                          region_name=settings.AWS_REGION_NAME,
                          aws_access_key_id=settings.AWS_ACCESS_KEY_ID,
                          aws_secret_access_key=settings.AWS_SECRET_ACCESS_KEY)

    try:
        # Provide the contents of the email.
        response = client.publish(
            TopicArn=settings.AWS_TOPIC_ARN,
            Message=instance.feedback,
            Subject=f'[ConnecTF] Feedback from: {instance.name}'
        )
    # Display an error if something goes wrong.
    except ClientError as e:
        logging.error(e.response['Error']['Message'])
    else:
        logging.info(f"Email sent! Message ID: {response['MessageId']}")


class FeedbackForm(ModelForm):
    token = CharField()

    def clean_token(self):
        token = self.cleaned_data['token']

        if token and settings.RECAPTCHA_SECRET:
            try:
                resp = requests.post(
                    'https://www.google.com/recaptcha/api/siteverify',
                    {
                        'secret': settings.RECAPTCHA_SECRET,
                        'response': token
                    })

                data = resp.json()

                if not (data['success'] and data['score'] >= 0.5):
                    raise ValidationError('you are a bot.')
            except requests.exceptions.RequestException as e:
                raise ValidationError("something wrong with recaptcha validation") from e

        return token

    class Meta:
        model = models.Feedback
        fields = ('name', 'feedback')


class FeedbackList(View):
    def post(self, request: HttpRequest):

        form = FeedbackForm(data=request.POST)
        if form.is_valid():
            if form.cleaned_data['token']:
                aws_send_email(form.instance)
            form.save()
            return HttpResponse(status=201)

        return HttpResponseBadRequest(form.errors)
