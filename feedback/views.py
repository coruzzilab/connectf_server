from django.forms.models import modelform_factory
from django.http import HttpRequest, HttpResponse, HttpResponseBadRequest
from django.views import View

from . import models

FeedbackForm = modelform_factory(models.Feedback, fields=('name', 'feedback'))


class FeedbackList(View):
    def post(self, request: HttpRequest):

        form = FeedbackForm(data=request.POST)
        if form.is_valid():
            form.save()
            return HttpResponse(status=201)

        return HttpResponseBadRequest()
