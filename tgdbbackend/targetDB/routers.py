from rest_framework.routers import DynamicDetailRoute, Route, SimpleRouter


########################################################################
class CustomReadOnlyRouter(SimpleRouter):
    """"""

    routes = [
        Route(
            url=r'^{prefix}{trailing_slash}$',
            mapping={'get': 'list'},
            name='{basename}-list',
            initkwargs={'suffix': 'List'}
        ),
        DynamicDetailRoute(
            url=r'^{prefix}/{methodnamehyphen}/{lookup}{trailing_slash}$',
            name='{basename}-{methodnamehyphen}',
            initkwargs={}
        ),
    ]
