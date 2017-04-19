from django.conf.urls import url
from . import views


app_name = 'upload'
urlpatterns = [
    url(r'^$', views.UploadView.as_view())
]
