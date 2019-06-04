"""connectf URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/1.11/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  url(r'^$', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  url(r'^$', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.conf.urls import url, include
    2. Add a URL to urlpatterns:  url(r'^blog/', include('blog.urls'))
"""
from django.conf import settings
from django.conf.urls.static import static
from django.urls import include, re_path
from django.views import defaults as default_views

urlpatterns = [
                  re_path(r"^api/", include(('targetdb.urls', 'targetdb'),
                                            namespace='targetdb')),
                  re_path(r'^api/', include(('querytgdb.urls', 'queryapp'),
                                            namespace='queryapp')),
                  re_path(r'^api/', include(('overview.urls', 'overview'),
                                            namespace='overview'))
              ] + static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)

if settings.DEBUG:
    # This allows the error pages to be debugged during development, just visit
    # these url in browser to see how these error pages look like.

    urlpatterns += [
        re_path(r'^400/$', default_views.bad_request,
                kwargs={'exception': Exception('Bad Request!')}),
        re_path(r'^403/$', default_views.permission_denied,
                kwargs={'exception': Exception('Permission Denied')}),
        re_path(r'^404/$', default_views.page_not_found,
                kwargs={'exception': Exception('Page not Found')}),
        re_path(r'^500/$', default_views.server_error)
    ]
