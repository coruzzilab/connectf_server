from django.urls import re_path

from . import views

urlpatterns = [
    re_path(r'^(?:(?P<request_id>\d{4}-\d{2}-\d{2}T\d{9}Z\d{3})/)?$', views.QueryView.as_view(), name="queryapp"),
    re_path(r'^network/(?P<request_id>\d{4}-\d{2}-\d{2}T\d{9}Z\d{3})/$', views.NetworkJSONView.as_view()),
    re_path(r'^stats/(?P<request_id>\d{4}-\d{2}-\d{2}T\d{9}Z\d{3})/$', views.StatsView.as_view()),
    re_path(r'^export/(?P<request_id>\d{4}-\d{2}-\d{2}T\d{9}Z\d{3})?\.zip$', views.FileExportView.as_view()),
    re_path(r'^list_enrichment/(?P<request_id>\d{4}-\d{2}-\d{2}T\d{9}Z\d{3})\.svg$',
            views.ListEnrichmentSVGView.as_view()),
    re_path(r'^list_enrichment/(?P<request_id>\d{4}-\d{2}-\d{2}T\d{9}Z\d{3})/$',
            views.ListEnrichmentTableView.as_view()),
    re_path(r'^list_enrichment/(?P<request_id>\d{4}-\d{2}-\d{2}T\d{9}Z\d{3})/legend/$',
            views.ListEnrichmentLegendView.as_view()),
    re_path(r'^motif_enrichment/cluster_info.csv$', views.MotifEnrichmentInfo.as_view()),
    re_path(r'^motif_enrichment/(?P<request_id>\d{4}-\d{2}-\d{2}T\d{9}Z\d{3})?/$',
            views.MotifEnrichmentJSONView.as_view()),
    re_path(r'^motif_enrichment/(?P<request_id>\d{4}-\d{2}-\d{2}T\d{9}Z\d{3})?/heatmap\.svg$',
            views.MotifEnrichmentHeatmapView.as_view()),
    re_path(r'^motif_enrichment/(?P<request_id>\d{4}-\d{2}-\d{2}T\d{9}Z\d{3})?/heatmap_table/$',
            views.MotifEnrichmentHeatmapTableView.as_view()),
    re_path(r'^analysis_enrichment/(?P<request_id>\d{4}-\d{2}-\d{2}T\d{9}Z\d{3})/$',
            views.AnalysisEnrichmentView.as_view()),
    re_path(r'^summary/(?P<request_id>\d{4}-\d{2}-\d{2}T\d{9}Z\d{3})/$', views.SummaryView.as_view()),
    re_path(r'^aupr/(?P<request_id>\d{4}-\d{2}-\d{2}T\d{9}Z\d{3})/$', views.NetworkAuprView.as_view())
]
