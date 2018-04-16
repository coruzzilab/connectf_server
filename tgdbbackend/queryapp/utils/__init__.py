import numpy as np
from django.core.serializers.json import DjangoJSONEncoder

__all__ = ('PandasJSONEncoder',)


class PandasJSONEncoder(DjangoJSONEncoder):
    def default(self, o):
        if np.isinf(o) or np.isnan(o):
            return None
        if isinstance(o, np.integer):
            return int(o)
        if isinstance(o, np.floating):
            return float(o)

        return super().default(o)
