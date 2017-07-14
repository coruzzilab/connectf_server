from django.conf.urls import url
from . import views

urlpatterns = [
    url(r'^$', views.FeedbackList.as_view()),
    url(r'^(?P<pk>[0-9]+)/$', views.FeedbackDetail.as_view())
]
