from django.urls import include, path, re_path

from targetdb import views

urlpatterns = [
    re_path(r'^feedback/', include('feedback.urls')),
    re_path(r'^edges/$', views.EdgeListView.as_view()),
    re_path(r'^lists/$', views.InterestingListsView.as_view()),
    re_path(r'^networks/$', views.InterestingNetworksView.as_view()),
    re_path(r'^tfs/$', views.TFView.as_view()),
    re_path(r'^key/$', views.KeyView.as_view()),
    re_path(r'^key/(?P<key>.+)/$', views.ValueView.as_view()),
    path('extra_field_names/', views.ExtraFieldNamesView.as_view())
]
