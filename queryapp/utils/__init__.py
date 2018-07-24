import gzip
import mimetypes
import pickle
from pathlib import Path
from typing import Union
from uuid import UUID

import numpy as np
import pandas as pd
from django.core.serializers.json import DjangoJSONEncoder


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


class CytoscapeJSONEncoder(DjangoJSONEncoder):
    def default(self, o):
        if isinstance(o, UUID):
            return str(o)
        return super().default(o)


def cache_result(obj, cache_name: Union[str, Path]):
    encoding = mimetypes.guess_type(cache_name)[1]
    if encoding == 'gzip':
        cache = gzip.open(cache_name, 'wb')
    else:
        cache = open(cache_name, 'wb')

    pickle.dump(obj, cache, protocol=pickle.HIGHEST_PROTOCOL)
    cache.close()
