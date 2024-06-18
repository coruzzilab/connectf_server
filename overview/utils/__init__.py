import re
from itertools import cycle, islice

import numpy as np
import pandas as pd
import pyparsing as pp
from scipy.interpolate import interp1d

name = pp.Word(pp.pyparsing_unicode.alphanums + '-_.:/\\') | pp.QuotedString('"', escChar='\\') | pp.QuotedString("'",
                                                                                                                  escChar='\\')
modname = name + pp.oneOf('= == !=') + name

na_vals = re.compile(r'^(?:N/?A?|NONE|NAN|\s*)$', flags=re.I)

COLORS = [
    '#a6cee3', '#1f78b4', '#b2df8a',
    '#33a02c', '#fb9a99', '#e31a1c',
    '#fdbf6f', '#ff7f00', '#cab2d6',
    '#6a3d9a', '#ffff99', '#b15928'
]

TYPE_COLORS = {
    "green": np.array(
        [[0, 68, 27], [27, 120, 55], [90, 174, 97], [166, 219, 160], [217, 240, 211]]),
    "purple": np.array(
        [[64, 0, 75], [118, 42, 131], [153, 112, 171], [194, 165, 207], [231, 212, 232]])
}

grfunc = interp1d(np.linspace(0, 1, num=5), TYPE_COLORS["green"][:, 0])
ggfunc = interp1d(np.linspace(0, 1, num=5), TYPE_COLORS["green"][:, 1])
gbfunc = interp1d(np.linspace(0, 1, num=5), TYPE_COLORS["green"][:, 2])

prfunc = interp1d(np.linspace(0, 1, num=5), TYPE_COLORS["purple"][:, 0])
pgfunc = interp1d(np.linspace(0, 1, num=5), TYPE_COLORS["purple"][:, 1])
pbfunc = interp1d(np.linspace(0, 1, num=5), TYPE_COLORS["purple"][:, 2])


def color_to_hex(c):
    return "#{1[0]:{0}}{1[1]:{0}}{1[2]:{0}}".format('02x', c)


class TypeColor:
    def __getitem__(self, item):
        if item == 'EXPRESSION':
            return self.get_expression_colors

        if item == 'BINDING':
            return self.get_binding_colors

        return self.get_default_colors

    def __call__(self, df: pd.DataFrame):
        type_count = df['type'].value_counts()
        colors = {str(idx): iter(self[idx](value)) for idx, value in type_count.items()}
        return pd.DataFrame(
            ((c, next(colors[str(t)])) for c, t in df.itertuples(index=False, name=None)),
            index=df.index,
            columns=['count', 'type'])

    @staticmethod
    def get_expression_colors(num: int):
        c = np.linspace(0, 1, num=num)

        return np.apply_along_axis(
            color_to_hex,
            1,
            np.concatenate([
                grfunc(c)[:, np.newaxis],
                ggfunc(c)[:, np.newaxis],
                gbfunc(c)[:, np.newaxis]
            ], axis=1).round().astype(int)
        )

    @staticmethod
    def get_binding_colors(num: int):
        c = np.linspace(0, 1, num=num)

        return np.apply_along_axis(
            color_to_hex,
            1,
            np.concatenate([
                prfunc(c)[:, np.newaxis],
                pgfunc(c)[:, np.newaxis],
                pbfunc(c)[:, np.newaxis]
            ], axis=1).round().astype(int)
        )

    @staticmethod
    def get_default_colors(num: int):
        return list(islice(cycle(COLORS), num))


TYPE_COLOR = TypeColor()
