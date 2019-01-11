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
        'color': '#FF9900',
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
            'INDUCED': '#4daf4a',
            'REPRESSED': '#e41a1c',
            'BOUND': '#377eb8'
        })

    def __getitem__(self, item):
        return super().__getitem__(simplify_edge(item))


COLOR = Color()

COLOR_SHAPE = defaultdict(default_shape, {
    'TXNFACTOR': {
        'color': '#00FF00',
        'shape': 'triangle'
    },
    'PROTEIN_CODING': {
        'color': '#AED6F1',
        'shape': 'roundrectangle'
    },
    'METABOLIC': {
        'color': '#D0ECE7',
        'shape': 'roundrectangle'
    },
    'MOLECULE': {
        'color': '#FF9900',
        'shape': 'roundrectangle'
    }
})
