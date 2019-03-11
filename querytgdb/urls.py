from django.urls import path, register_converter

import sungear.views
from . import views


class FloatConverter:
    regex = r"(?:[-+]?\d*\.?\d+)(?:[eE]([-+]?\d+))?"

    def to_python(self, value):
        return float(value)

    def to_url(self, value):
        return str(value)


register_converter(FloatConverter, 'float')

urlpatterns = [
    path('', views.QueryView.as_view(), name="queryapp"),
    path('<uuid:request_id>/', views.QueryView.as_view()),
    path('network/<uuid:request_id>/', views.NetworkJSONView.as_view()),
    path('stats/<uuid:request_id>/', views.StatsView.as_view()),
    path('export/<uuid:request_id>.zip', views.FileExportView.as_view()),
    path('list_enrichment/<uuid:request_id>.svg',
         views.ListEnrichmentSVGView.as_view()),
    path('list_enrichment/<uuid:request_id>/',
         views.ListEnrichmentTableView.as_view()),
    path('list_enrichment/<uuid:request_id>/legend/',
         views.ListEnrichmentLegendView.as_view()),
    path('motif_enrichment/cluster_info.csv', views.MotifEnrichmentInfo.as_view()),
    path('motif_enrichment/<uuid:request_id>/',
         views.MotifEnrichmentJSONView.as_view()),
    path('motif_enrichment/<uuid:request_id>/heatmap.svg',
         views.MotifEnrichmentHeatmapView.as_view()),
    path('motif_enrichment/<uuid:request_id>/heatmap_table/',
         views.MotifEnrichmentHeatmapTableView.as_view()),
    path('analysis_enrichment/<uuid:request_id>/',
         views.AnalysisEnrichmentView.as_view()),
    path('summary/<uuid:request_id>/', views.SummaryView.as_view()),
    path('aupr/<uuid:request_id>/', views.NetworkAuprView.as_view()),
    path('aupr/<uuid:request_id>/pruned/<float:cutoff>/', views.NetworkPrunedView.as_view())
]
