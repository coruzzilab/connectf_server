from django.conf.urls import include, url

from tgdbbackend.targetDB import views
from .routers import CustomReadOnlyRouter

# router = routers.DefaultRouter()
router = CustomReadOnlyRouter()
router.register(r'metas', views.MetaValueDistinctViewSet)
router.register(r'tfs', views.TFValueDistinctViewSet)
router.register(r'edges', views.EdgesValueDistinctViewSet)

urlpatterns = [
    url(r'^feedback/', include('feedback.urls'))
] + router.urls
