import base64
import io
import logging
import math
import pkgutil
import sys
from collections import UserDict
from concurrent.futures import Future, ThreadPoolExecutor
from functools import wraps
from operator import methodcaller
from typing import Any, Callable, Dict, Iterable, Optional, Sized, TypeVar
from uuid import UUID

import numpy as np
import pandas as pd
from django.core.serializers.json import DjangoJSONEncoder
from django.db import DatabaseError
from django.db.models import QuerySet
from django.http import FileResponse
from fontTools.ttLib import TTFont
from lxml import etree

from querytgdb.models import Analysis, AnalysisData, Annotation

logger = logging.getLogger(__name__)

T = TypeVar('T')


class PandasJSONEncoder(DjangoJSONEncoder):
    def default(self, o):
        if isinstance(o, np.number):
            if np.isinf(o):
                return math.inf
            if np.isnan(o):
                return math.nan
            return o.item()
        if isinstance(o, (pd.Index, pd.Series, np.ndarray)):
            return o.tolist()

        return super().default(o)


class NetworkJSONEncoder(DjangoJSONEncoder):
    def default(self, o):
        if isinstance(o, UUID):
            return str(o)
        return super().default(o)


class GzipFileResponse(FileResponse):
    """
    Handles gzip files correctly

    Does not compress files
    """

    def __init__(self, *args, as_attachment=False, filename='', **kwargs):
        super().__init__(*args, as_attachment=as_attachment, filename=filename, **kwargs)
        self['Content-Encoding'] = 'gzip'


class CaselessDict(UserDict):
    def __init__(self, dict_, **kwargs):
        super().__init__({k.lower(): v for k, v in dict_.items()}, **kwargs)

    def __getitem__(self, item):
        return super().__getitem__(item.lower())

    def __setitem__(self, key, value):
        return super().__setitem__(key.lower(), value)


def convert_float(s) -> Optional[float]:
    try:
        return float(s)
    except (ValueError, TypeError):
        return None


def metadata_to_dict(df: pd.DataFrame) -> Dict[str, Any]:
    df = df.reset_index()
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


with TTFont(io.BytesIO(pkgutil.get_data('matplotlib', 'mpl-data/fonts/ttf/DejaVuSans.ttf')),
            recalcBBoxes=False,
            recalcTimestamp=False) as font, \
        io.BytesIO() as font_buff:
    font.flavor = "woff"
    font.save(font_buff, reorderTables=False)
    FONT_STR = '@font-face {{font-family: "DejaVu Sans"; src: local("DejaVu Sans"), local("DejaVuSans"), ' \
               'url({}) format("woff");}}'.format('data:font/woff;base64,' +
                                                  base64.standard_b64encode(font_buff.getvalue()).decode())


def svg_font_adder(buff: io.BytesIO) -> io.BytesIO:
    """
    Adds font file as a base64 encoded string to an SVG

    The very definition of overkill
    :param buff:
    :return:
    """
    tree = etree.parse(buff)
    style = tree.find('./{http://www.w3.org/2000/svg}defs/{http://www.w3.org/2000/svg}style')

    style.text += FONT_STR

    buff.truncate(0)
    tree.write(buff)
    buff.seek(0)

    return buff


def clear_data(df: pd.DataFrame, drop: bool = True) -> pd.DataFrame:
    """
    Remove unneeded columns when doing later calculations
    :param df:
    :param drop:
    :return:
    """
    df = df.loc[:, (slice(None), slice(None), ['EDGE', 'Log2FC'])]
    if drop:
        df.columns = df.columns.droplevel(2)

    return df


def data_to_edges(df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert EDGE and Log2FC to respective edge_type
    :param df:
    :return:
    """
    df = df.loc[:, (slice(None), slice(None), ['EDGE', 'Log2FC'])]

    edge_types = dict(AnalysisData.objects.filter(
        key__name='EDGE_TYPE',
        analysis_id__in=df.columns.get_level_values(1)
    ).values_list('analysis_id', 'value'))

    def set_edge_name(s):
        try:
            edge_type = edge_types[s.name[1]]
        except KeyError:
            edge_type = 'edge'

        c = pd.Series(index=s.index, dtype=str)

        if s.name[2] == 'Log2FC':
            c = c.mask(s >= 0, edge_type + ':INDUCED')
            c = c.mask(s < 0, edge_type + ':REPRESSED')

            return c
        elif (s == "*").any():
            return s.mask(s.notna(), edge_type + ':EXPRESSION')
        else:
            return s.mask(s.notna(), edge_type)

    df = df.apply(set_edge_name)

    df.columns = df.columns.droplevel(2)

    return df


def get_size(func: Callable[..., Sized]) -> Callable[..., Sized]:
    """
    Get size of function out put

    :param func:
    :return:
    """

    @wraps(func)
    def wrapper(*args, **kwargs):
        result = func(*args, **kwargs)

        logger.info(f"{func.__name__} len: {len(result)} size: {sys.getsizeof(result)}")

        return result

    return wrapper


class AsyncDataLoader:
    def __init__(self):
        self.data = {}
        self.pool = ThreadPoolExecutor()

    def __setitem__(self, key, value):
        if callable(value):
            self.data[key] = self.pool.submit(value)
        else:
            self.data[key] = value

    def __getitem__(self, item):
        if isinstance(self.data[item], Future):
            result = self.data[item].result()
            self.data[item] = result
            return result

        return self.data[item]

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def close(self, *args, **kwargs):
        self.pool.shutdown(*args, **kwargs)


def skip_for_management(func):
    """
    Return a noop when not run in WSGI
    :param func:
    :return:
    """
    if 'django.core.wsgi' in sys.modules:
        return func

    @wraps(func)
    def f(*args, **kwargs):
        return None

    return f


def get_annotations():
    try:
        anno = pd.DataFrame(
            Annotation.objects.values_list(
                'gene_id', 'fullname', 'gene_family', 'gene_type', 'name', 'id').iterator(),
            columns=['TARGET', 'Full Name', 'Gene Family', 'Type', 'Name', 'id'])
        anno = anno.set_index('TARGET')
    except DatabaseError:
        anno = pd.DataFrame(columns=['Full Name', 'Gene Family', 'Type', 'Name', 'id'])
        anno.index.name = 'TARGET'

    return anno


async_loader = AsyncDataLoader()
async_loader['annotations'] = get_annotations


def check_annotations(genes):
    return set(map(methodcaller('upper'), genes)) - set(async_loader['annotations'].index.str.upper())


def get_metadata(analyses, fields: Optional[Iterable[str]] = None) -> pd.DataFrame:
    if not isinstance(analyses, QuerySet):
        analyses = Analysis.objects.filter(pk__in=analyses)

    opts = {}

    if fields:
        opts['key__name__in'] = fields

    analysis_data = AnalysisData.objects.filter(analysis__in=analyses, **opts).prefetch_related(
        'key')
    metadata = pd.DataFrame(
        analysis_data.values_list('analysis_id', 'key__name', 'value').iterator(chunk_size=2000),
        columns=['id', 'key', 'value'])
    metadata = metadata.set_index(['id', 'key'])['value'].unstack().fillna('None')

    genes = pd.DataFrame(analyses.values_list('id', 'tf__gene_id', 'tf__name').iterator(),
                         columns=['id', 'gene_id', 'gene_name'])
    genes = genes.set_index('id')
    metadata = genes.merge(metadata, left_index=True, right_index=True)

    metadata.insert(0, 'analysis_id', metadata.index.astype(str))

    if fields:
        metadata = metadata.reindex(columns=fields)

    metadata = metadata.fillna('')

    return metadata
