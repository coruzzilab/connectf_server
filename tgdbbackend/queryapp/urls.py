from django.conf.urls import url

from tgdbbackend.queryapp import views

urlpatterns = [
    url(r'^$', views.HandleQueryView.as_view(), name="queryapp"),
]
