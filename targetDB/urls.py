from django.urls import include, re_path

from targetDB import views

# from .routers import CustomReadOnlyRouter
#
# router = CustomReadOnlyRouter()
# router.register(r'metas', views.MetaValueDistinctViewSet)
# router.register(r'tfs', views.TFValueDistinctViewSet)
# router.register(r'edges', views.EdgesValueDistinctViewSet)

urlpatterns = [
    re_path(r'^feedback/', include('feedback.urls')),
    re_path(r'^lists/$', views.InterestingListsView.as_view()),
    re_path(r'^tfs/$', views.TFView.as_view()),
    re_path(r'^key/$', views.KeyView.as_view()),
    re_path(r'^key/(?P<key>.+)/$', views.ValueView.as_view())
]
