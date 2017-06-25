from django.conf.urls import url

from tgdbbackend.queryapp import views

urlpatterns = [
    url(r'^$', views.HandleQueryView.as_view(), name="queryapp"),
    url(r'^cytoscape/(?P<request_id>\d{4}-\d{2}-\d{2}T\d{9}Z\d{1,3})/(?P<name>[a-egimnorstvwy1-3_]+)/$',
        views.CytoscapeJSONView.as_view()),
    url(r'^excel/(?P<request_id>\d{4}-\d{2}-\d{2}T\d{9}Z\d{1,3})\.zip', views.ExcelDownloadView.as_view()),
    url(r'^heatmap/(?P<request_id>\d{4}-\d{2}-\d{2}T\d{9}Z\d{1,3})\.png', views.HeatMapPNGView.as_view())
]
