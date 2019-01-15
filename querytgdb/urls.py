from django.urls import re_path

from . import views

urlpatterns = [
    re_path(r'^(?:(?P<request_id>[a-f0-9\-]+)/)?$', views.QueryView.as_view(), name="queryapp"),
    re_path(r'^network/(?P<request_id>[a-f0-9\-]+)/$', views.NetworkJSONView.as_view()),
    re_path(r'^stats/(?P<request_id>[a-f0-9\-]+)/$', views.StatsView.as_view()),
    re_path(r'^export/(?P<request_id>[a-f0-9\-]+)?\.zip$', views.FileExportView.as_view()),
    re_path(r'^list_enrichment/(?P<request_id>[a-f0-9\-]+)\.svg$',
            views.ListEnrichmentSVGView.as_view()),
    re_path(r'^list_enrichment/(?P<request_id>[a-f0-9\-]+)/$',
            views.ListEnrichmentTableView.as_view()),
    re_path(r'^list_enrichment/(?P<request_id>[a-f0-9\-]+)/legend/$',
            views.ListEnrichmentLegendView.as_view()),
    re_path(r'^motif_enrichment/cluster_info.csv$', views.MotifEnrichmentInfo.as_view()),
    re_path(r'^motif_enrichment/(?P<request_id>[a-f0-9\-]+)?/$',
            views.MotifEnrichmentJSONView.as_view()),
    re_path(r'^motif_enrichment/(?P<request_id>[a-f0-9\-]+)?/heatmap\.svg$',
            views.MotifEnrichmentHeatmapView.as_view()),
    re_path(r'^motif_enrichment/(?P<request_id>[a-f0-9\-]+)?/heatmap_table/$',
            views.MotifEnrichmentHeatmapTableView.as_view()),
    re_path(r'^analysis_enrichment/(?P<request_id>[a-f0-9\-]+)/$',
            views.AnalysisEnrichmentView.as_view()),
    re_path(r'^summary/(?P<request_id>[a-f0-9\-]+)/$', views.SummaryView.as_view()),
    re_path(r'^aupr/(?P<request_id>[a-f0-9\-]+)/$', views.NetworkAuprView.as_view())
]
