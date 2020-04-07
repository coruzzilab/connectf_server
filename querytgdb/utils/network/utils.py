from collections import UserDict, defaultdict
from functools import lru_cache

__all__ = ['COLOR', 'COLOR_SHAPE']


def default_shape(*args):
    """
    Default shape for genes of unknown type.

    :param args:
    :return:
    """
    return {
        'color': '#CCCCFF',
        'shape': 'roundrectangle'
    }


@lru_cache()
def simplify_edge(edge: str) -> str:
    if edge.endswith("INDUCED"):
        return "INDUCED"
    elif edge.endswith("REPRESSED"):
        return "REPRESSED"
    return "BOUND"


class Color(UserDict):
    def __init__(self):
        super().__init__({
            'INDUCED': {'color': '#4daf4a', 'shape': 'triangle'},
            'REPRESSED': {'color': '#e41a1c', 'shape': 'tee'},
            'BOUND': {'color': '#377eb8', 'shape': 'triangle'}
        })

    def __getitem__(self, item):
        return super().__getitem__(simplify_edge(item))


COLOR = Color()

COLOR_SHAPE = defaultdict(default_shape, {
    'TXNFACTOR': {
        'color': '#00CC00',
        'shape': 'triangle'
    },
    'PROTEIN_CODING': {
        'color': '#CCCCFF',
        'shape': 'roundrectangle'
    },
    'METABOLIC': {
        'color': '#3399FF',
        'shape': 'roundrectangle'
    },
    'MOLECULE': {
        'color': '#FFFF00',
        'shape': 'roundrectangle'
    }
})
