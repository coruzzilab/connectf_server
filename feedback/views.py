from rest_framework import generics

from . import models, serializers


class FeedbackList(generics.ListCreateAPIView):
    queryset = models.Feedback.objects.all()
    serializer_class = serializers.FeedbackSerializer


class FeedbackDetail(generics.RetrieveDestroyAPIView):
    queryset = models.Feedback.objects.all()
    serializer_class = serializers.FeedbackSerializer
