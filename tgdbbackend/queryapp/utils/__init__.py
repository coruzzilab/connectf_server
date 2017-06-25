import numpy as np
from django.core.serializers.json import DjangoJSONEncoder

__all__ = ('PandasJSONEncoder',)


class PandasJSONEncoder(DjangoJSONEncoder):
    def default(self, o):
        if isinstance(o, np.integer):
            return int(o)
        else:
            return super().default(o)
