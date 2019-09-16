import base64
import io
import logging
import math
import pkgutil
import re
import sys
from functools import wraps
from operator import methodcaller
from threading import Lock, Thread
from typing import Any, Callable, Dict, Iterable, Optional, Set, Sized, Tuple, TypeVar
from uuid import UUID

import numpy as np
import pandas as pd
from django.core.serializers.json import DjangoJSONEncoder
from django.db import DatabaseError
from django.http import FileResponse
from fontTools.ttLib import TTFont
from lxml import etree

from querytgdb.models import AnalysisData, Annotation

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


NAME_REGEX = re.compile(
    r"^([^\"]+)(?:\s*\"([^\"]+)\")?\s*([a-f0-9]{8}-[a-f0-9]{4}-4[a-f0-9]{3}-[89aAbB][a-f0-9]{3}-[a-f0-9]{12})?",
    flags=re.I)


def split_name(name: str) -> Tuple[str, str, str]:
    m = NAME_REGEX.match(name)

    if m:
        name, criterion, uid = map(methodcaller('strip'), m.groups(''))
    else:
        criterion = uid = ''

    return name, criterion, uid


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


class Annotations:
    def __init__(self):
        self._anno: pd.DataFrame = None

        self.lock = Lock()
        self.task = Thread(target=self.get_annotations)
        self.task.start()

    def get_annotations(self):
        """
        Loads annotations from the database into memory.
        """
        try:
            anno = pd.DataFrame(
                Annotation.objects.values_list(
                    'gene_id', 'fullname', 'gene_family', 'gene_type', 'name', 'id').iterator(),
                columns=['TARGET', 'Full Name', 'Gene Family', 'Type', 'Name', 'id'])
            anno = anno.set_index('TARGET')
        except DatabaseError:
            anno = pd.DataFrame(columns=['Full Name', 'Gene Family', 'Type', 'Name', 'id'])
            anno.index.name = 'TARGET'

        with self.lock:
            self._anno = anno

    def __call__(self):
        return self.annotation

    @property
    def annotation(self) -> pd.DataFrame:
        with self.lock:
            if self._anno is None:
                self.task.join()

        return self._anno

    @property
    def genes(self) -> Set:
        return set(self.annotation.index)

    @property
    def genes_upper(self) -> Set:
        return set(annotations.annotation.index.str.upper())


annotations = Annotations()


def check_annotations(genes):
    return set(map(methodcaller('upper'), genes)) - annotations.genes_upper


def get_metadata(analyses, fields: Iterable[str] = None) -> pd.DataFrame:
    opts = {}

    if fields is not None:
        opts['key__name__in'] = fields

    analysis_data = AnalysisData.objects.filter(analysis__in=analyses, **opts).prefetch_related(
        'key')
    metadata = pd.DataFrame(
        analysis_data.values_list('analysis_id', 'key__name', 'value').iterator(),
        columns=['id', 'key', 'value'])
    metadata = metadata.set_index(['id', 'key'])['value'].unstack().fillna('None')

    genes = pd.DataFrame(analyses.values_list('id', 'tf__gene_id', 'tf__name').iterator(),
                         columns=['id', 'gene_id', 'gene_name'])
    genes = genes.set_index('id')
    metadata = genes.merge(metadata, left_index=True, right_index=True)

    metadata.insert(0, 'analysis_id', metadata.index.astype(str))

    if fields:
        metadata = metadata.loc[:, fields]

    return metadata
