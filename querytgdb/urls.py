from django.urls import path

from . import views

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
    path('aupr/<uuid:request_id>/pruned/(P<cutoff>[+\-0-9eE.]+)/', views.NetworkPrunedView.as_view())
]
