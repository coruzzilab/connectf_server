from rest_framework.routers import DynamicDetailRoute, Route, SimpleRouter


########################################################################
class CustomReadOnlyRouter(SimpleRouter):
    """"""

    routes = [
        Route(
            url=r'^{prefix}$',
            mapping={'get': 'list'},
            name='{basename}-list',
            initkwargs={'suffix': 'List'}
        ),
        DynamicDetailRoute(
            url=r'^{prefix}/{methodnamehyphen}/{lookup}$',
            name='{basename}-{methodnamehyphen}',
            initkwargs={}
        ),
    ]
