import numpy as np
from django.core.serializers.json import DjangoJSONEncoder
import pandas as pd

__all__ = ('PandasJSONEncoder',)


class PandasJSONEncoder(DjangoJSONEncoder):
    def default(self, o):
        if isinstance(o, np.number):
            if np.isinf(o) or np.isnan(o):
                return None
            if isinstance(o, np.integer):
                return int(o)
            if isinstance(o, np.floating):
                return float(o)
        if isinstance(o, pd.Index):
            return o.tolist()

        return super().default(o)
