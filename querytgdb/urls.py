from django.urls import path, register_converter

import sungear_app.views
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
    path('ids/<uuid:request_id>/', views.EditQueryView.as_view()),
    path('network/<uuid:request_id>/', views.NetworkJSONView.as_view()),
    path('network/<uuid:request_id>.sif', views.NetworkSifView.as_view()),
    path('stats/<uuid:request_id>/', views.StatsView.as_view()),
    path('export/<uuid:request_id>.xlsx', views.ExcelExportView.as_view()),
    path('export/<uuid:request_id>.csv', views.CsvExportView.as_view()),
    path('list_enrichment/<uuid:request_id>.svg',
         views.ListEnrichmentHeatmapView.as_view()),
    path('list_enrichment/<uuid:request_id>/',
         views.ListEnrichmentTableView.as_view()),
    path('list_enrichment/<uuid:request_id>/legend/',
         views.ListEnrichmentLegendView.as_view()),
    path('motif_enrichment/cluster_info.csv', views.MotifEnrichmentInfoView.as_view()),
    path('motif_enrichment/regions/', views.MotifEnrichmentRegionsView.as_view()),
    path('motif_enrichment/motifs/', views.MotifEnrichmentMotifsView.as_view()),
    path('motif_enrichment/additional/motifs/', views.MotifEnrichmentAdditionalMotifsView.as_view()),
    path('motif_enrichment/<uuid:request_id>/',
         views.MotifEnrichmentJSONView.as_view()),
    path('motif_enrichment/additional/<uuid:request_id>/',
         views.AdditionalMotifEnrichmentJSONView.as_view()),
    path('motif_enrichment/<uuid:request_id>/heatmap.svg',
         views.MotifEnrichmentHeatmapView.as_view()),
    path('motif_enrichment/<uuid:request_id>/heatmap_table/',
         views.MotifEnrichmentHeatmapTableView.as_view()),
    path('analysis_enrichment/<uuid:request_id>/',
         views.AnalysisEnrichmentView.as_view()),
    path('analysis_enrichment/<uuid:request_id>.csv', views.AnalysisEnrichmentCsvView.as_view()),
    path('summary/<uuid:request_id>/', views.SummaryView.as_view()),
    path('aupr/<uuid:request_id>/', views.NetworkAuprView.as_view()),
    path('aupr/<uuid:request_id>/pruned/<float:cutoff>/', views.NetworkPrunedView.as_view()),
    path('sungear/<uuid:request_id>/', sungear_app.views.SungearView.as_view()),
    path('list_download/<str:list_name>/', views.ListDownloadView.as_view())
]
