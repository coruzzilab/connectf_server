from __future__ import absolute_import, unicode_literals

from django.conf.urls import url

from tgdbbackend.queryapp import views

urlpatterns = [
    url(r'^$', views.HandleQueryView.as_view(), name="queryapp"),
]
