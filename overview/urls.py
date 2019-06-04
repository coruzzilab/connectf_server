from django.urls import re_path

from . import views

urlpatterns = [
    re_path(r'^overview/$', views.OverviewView.as_view()),
    re_path(r'^overview/autocomplete/$', views.OverviewAutocompleteView.as_view())
]
