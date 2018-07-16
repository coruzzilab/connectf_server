from rest_framework.routers import DynamicDetailRoute, Route, SimpleRouter, DynamicListRoute


########################################################################
class CustomReadOnlyRouter(SimpleRouter):
    """"""

    routes = [
        Route(
            url=r'^{prefix}{trailing_slash}$',
            mapping={
                'get': 'list',
                'post': 'create'
            },
            name='{basename}-list',
            initkwargs={'suffix': 'List'},
            detail=False
        ),
        DynamicListRoute(
            url=r'^{prefix}/{methodname}{trailing_slash}$',
            name='{basename}-{methodnamehyphen}',
            initkwargs={},
        ),
        Route(
            url=r'^{prefix}{trailing_slash}$',
            mapping={'get': 'list'},
            name='{basename}-list',
            initkwargs={'suffix': 'List'},
            detail=False
        ),
        DynamicDetailRoute(
            url=r'^{prefix}/{methodname}/{lookup}{trailing_slash}$',
            name='{basename}-{methodnamehyphen}',
            initkwargs={}
        ),
    ]
