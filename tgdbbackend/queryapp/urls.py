from django.conf.urls import url

from tgdbbackend.queryapp import views

urlpatterns = [
    url(r'^$', views.HandleQueryView.as_view(), name="queryapp"),
    url(r'^cytoscape/(?P<request_id>\d{4}-\d{2}-\d{2}T\d{9}Z\d{3})/(?P<name>[a-egimnorstvwy1-3_]+)/$',
        views.CytoscapeJSONView.as_view()),
    url(r'^excel/(?P<request_id>\d{4}-\d{2}-\d{2}T\d{9}Z\d{3})\.zip$', views.ExcelDownloadView.as_view()),
    url(r'^heatmap/(?P<request_id>\d{4}-\d{2}-\d{2}T\d{9}Z\d{3})\.svg$', views.HeatMapPNGView.as_view()),
    url(r'^heatmap/(?P<request_id>\d{4}-\d{2}-\d{2}T\d{9}Z\d{3})/$', views.HeatMapTableView.as_view()),
    url(r'^motif_enrichment/(?P<request_id>\d{4}-\d{2}-\d{2}T\d{9}Z\d{3})/$', views.MotifEnrichmentView.as_view()),
    url(r'^motif_enrichment/(?P<request_id>\d{4}-\d{2}-\d{2}T\d{9}Z\d{3})/heatmap\.svg$',
        views.MotifEnrichmentHeatmapView.as_view()),
    url(r'^motif_enrichment/(?P<request_id>.*)/heatmap\.svg$', views.MotifEnrichmentErrorView.as_view())
]
