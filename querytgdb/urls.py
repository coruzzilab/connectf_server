from django.conf.urls import url

from . import views

urlpatterns = [
    url(r'^$', views.QueryView.as_view(), name="queryapp"),
    url(r'^(?P<request_id>\d{4}-\d{2}-\d{2}T\d{9}Z\d{3})/$', views.QueryIdView.as_view()),
    url(r'^cytoscape/(?P<request_id>\d{4}-\d{2}-\d{2}T\d{9}Z\d{3})/$', views.CytoscapeJSONView.as_view()),
    url(r'^stats/(?P<request_id>\d{4}-\d{2}-\d{2}T\d{9}Z\d{3})/$', views.StatsView.as_view()),
    url(r'^export/(?P<request_id>\d{4}-\d{2}-\d{2}T\d{9}Z\d{3})?\.zip$', views.FileExportView.as_view()),
    url(r'^list_enrichment/(?P<request_id>\d{4}-\d{2}-\d{2}T\d{9}Z\d{3})\.svg$', views.ListEnrichmentSVGView.as_view()),
    url(r'^list_enrichment/(?P<request_id>\d{4}-\d{2}-\d{2}T\d{9}Z\d{3})/$', views.ListEnrichmentTableView.as_view()),
    url(r'^list_enrichment/(?P<request_id>\d{4}-\d{2}-\d{2}T\d{9}Z\d{3})/legend/$',
        views.ListEnrichmentLegendView.as_view()),
    url(r'^motif_enrichment/(?P<request_id>\d{4}-\d{2}-\d{2}T\d{9}Z\d{3})?/$', views.MotifEnrichmentJSONView.as_view()),
    url(r'^motif_enrichment/(?P<request_id>\d{4}-\d{2}-\d{2}T\d{9}Z\d{3})?/heatmap\.svg$',
        views.MotifEnrichmentHeatmapView.as_view()),
    url(r'^motif_enrichment/(?P<request_id>\d{4}-\d{2}-\d{2}T\d{9}Z\d{3})?/heatmap_table/$',
        views.MotifEnrichmentHeatmapTableView.as_view()),
    url(r'^analysis_enrichment/(?P<request_id>\d{4}-\d{2}-\d{2}T\d{9}Z\d{3})/$', views.AnalysisEnrichmentView.as_view())
]
