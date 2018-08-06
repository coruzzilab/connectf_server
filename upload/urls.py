from django.conf.urls import url
from . import views


app_name = 'upload'
urlpatterns = [
    # url(r'^$', views.UploadExperimentView.as_view()),
    url(r'^analysis/$', views.UploadAnalysisView.as_view())
]
