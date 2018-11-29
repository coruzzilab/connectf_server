import base64
import gzip
import math
import mimetypes
import os.path as path
import pickle
import re
from contextlib import closing
from operator import methodcaller
from pathlib import Path
from typing import Any, Callable, Dict, Optional, Tuple, Union
from uuid import UUID

import numpy as np
import pandas as pd
from django.conf import settings
from django.core.exceptions import ObjectDoesNotExist
from django.core.serializers.json import DjangoJSONEncoder
from django.db.models import QuerySet
from lxml import etree

from ..models import Analysis


class PandasJSONEncoder(DjangoJSONEncoder):
    def default(self, o):
        if isinstance(o, np.number):
            if np.isinf(o):
                return math.inf
            if np.isnan(o):
                return None
            if isinstance(o, np.integer):
                return int(o)
            if isinstance(o, np.floating):
                return float(o)
        if isinstance(o, (pd.Index, pd.Series, np.ndarray)):
            return o.tolist()

        return super().default(o)


class CytoscapeJSONEncoder(DjangoJSONEncoder):
    def default(self, o):
        if isinstance(o, UUID):
            return str(o)
        return super().default(o)


def cache_result(obj, cache_name: Union[str, Path], *args, **kwargs) -> Any:
    """
    Caches object as gzipped picked
    :param obj:
    :param cache_name:
    :return:
    """
    encoding = mimetypes.guess_type(cache_name)[1]
    if encoding == 'gzip':
        cache = gzip.open(cache_name, 'wb', *args, **kwargs)
    else:
        cache = open(cache_name, 'wb')

    with closing(cache) as c:
        pickle.dump(obj, c, protocol=pickle.HIGHEST_PROTOCOL)


def read_cached_result(cache_name: Union[str, Path]) -> Any:
    """
    Read cached pickle
    :param cache_name:
    :return:
    """
    encoding = mimetypes.guess_type(cache_name)[1]
    if encoding == 'gzip':
        cache = gzip.open(cache_name, 'rb')
    else:
        cache = open(cache_name, 'rb')

    with closing(cache) as c:
        return pickle.load(c)


def cache_view(func: Callable, cache_path) -> Any:
    """
    read results from cache if possible
    :param func:
    :param cache_path:
    :return:
    """
    try:
        result = read_cached_result(cache_path)
    except FileNotFoundError:
        result = func()
        cache_result(result, cache_path)

    return result


def convert_float(s) -> Optional[float]:
    try:
        return float(s)
    except (ValueError, TypeError):
        return None


def metadata_to_dict(df: pd.DataFrame) -> Dict[str, Any]:
    df = df.rename(columns=lambda x: f"ID: {x}").reset_index()
    meta_columns = [{"id": column, "name": column, "field": column} for column in df.columns]

    df = df.where(df.notna(), None)

    return {'columns': meta_columns,
            'data': df.to_dict(orient='index')}


def get_exp_type(metadata: pd.Series) -> str:
    if metadata['EXPERIMENT_TYPE'] == 'EXPRESSION':
        return metadata['EXPRESSION_TYPE']
    elif metadata['EXPERIMENT_TYPE'] == 'BINDING':
        return metadata['BINDING_TYPE']
    else:
        raise ValueError('Invalid EXPERIMENT_TYPE')


def column_string(n: int) -> str:
    s = ""
    while n > 0:
        n, remainder = divmod(n - 1, 26)
        s = chr(65 + remainder) + s
    return s


def svg_font_adder(buff):
    """
    Adds font file as a base64 encoded string to an SVG

    The very definition of overkill
    :param buff:
    :return:
    """
    tree = etree.parse(buff)
    style = tree.find('./{http://www.w3.org/2000/svg}defs/{http://www.w3.org/2000/svg}style')

    with open(path.join(settings.BASE_DIR, 'tgdbbackend/static/fonts/DejaVuSans.woff'), 'rb') as font_file:
        font_str = base64.b64encode(font_file.read()).decode()

    style.text += '@font-face {{font-family: "DejaVu Sans"; src: local("DejaVu Sans"), local("DejaVuSans"), ' \
                  'url({}) format("woff");}}'.format('data:font/woff;base64,' + font_str)

    buff.truncate(0)
    tree.write(buff)

    return buff


NAME_REGEX = re.compile(r"^([^\"]+)(?:\s*\"([^\"]+)\")?")


def split_name(name: str) -> Tuple[str, str]:
    m = NAME_REGEX.match(name)

    if m:
        name, criterion = map(methodcaller('strip'), m.groups(''))
    else:
        criterion = ''

    return name, criterion


def clear_data(df: pd.DataFrame) -> pd.DataFrame:
    """
    Remove unneeded columns when doing later calculations
    :param df:
    :return:
    """
    df = df.loc[:, (slice(None), slice(None), ['EDGE', 'Log2FC'])]
    df.columns = df.columns.droplevel(2)

    return df


def data_to_edges(df: pd.DataFrame, analyses: Optional[QuerySet] = None, drop: bool = True) -> pd.DataFrame:
    """
    Convert EDGE and Log2FC to respective edge_type
    :param df:
    :param analyses:
    :param drop:
    :return:
    """
    df = df.loc[:, (slice(None), slice(None), ['EDGE', 'Log2FC'])]

    if analyses is None:
        analyses = Analysis.objects.filter(
            pk__in=df.columns.get_level_values(1)
        ).prefetch_related('analysisdata_set', 'analysisdata_set__key')

    for name, column in df.iteritems():
        try:
            edge_type = analyses.get(pk=name[1]).analysisdata_set.get(key__name='EDGE_TYPE').value
        except ObjectDoesNotExist:
            edge_type = 'edge'

        if name[2] == 'Log2FC':
            c = column.mask(column >= 0, edge_type + ':INDUCED')
            c = c.mask(column < 0, edge_type + ':REPRESSED')

            df.loc[:, name] = c
        else:
            df.loc[:, name] = column.mask(column.notna(), edge_type)

    if drop:
        df.columns = df.columns.droplevel(2)

    return df
