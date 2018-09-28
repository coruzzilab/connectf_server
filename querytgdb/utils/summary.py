import numpy as np
import pandas as pd
from collections import OrderedDict


def get_summary(cache_path):
    df = pd.read_pickle(cache_path)

    df = df.loc[:, (slice(None), slice(None), slice(None), 'EDGE')]
    df.columns = df.columns.droplevel(3)

    chart = df.apply(lambda x: x.value_counts()).fillna(0).astype(np.int_).to_dict(into=OrderedDict)

    chart = OrderedDict((','.join(key), value) for key, value in chart.items())

    return {
        'chart': chart
    }
