from rest_framework import generics

from . import models, serializers


class FeedbackList(generics.CreateAPIView):
    queryset = models.Feedback.objects.all()
    serializer_class = serializers.FeedbackSerializer
