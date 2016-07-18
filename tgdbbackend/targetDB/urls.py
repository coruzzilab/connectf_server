from __future__ import absolute_import, unicode_literals

from django.conf.urls import include, url

from rest_framework import routers;
from tgdbbackend.targetDB import views;
from rest_framework.urlpatterns import format_suffix_patterns
from .routers import CustomReadOnlyRouter

#router = routers.DefaultRouter();
router = CustomReadOnlyRouter();
router.register(r'metas',views.MetaValueDistinctViewSet);

urlpatterns = [
    url(r'^',include(router.urls)),
]