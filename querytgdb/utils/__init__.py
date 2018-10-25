import base64
import gzip
import math
import mimetypes
import os.path as path
import pickle
import re
import sys
from contextlib import closing
from pathlib import Path
from typing import Any, Dict, Optional, Tuple, Union
from uuid import UUID

import numpy as np
import pandas as pd
from django.conf import settings
from django.core.exceptions import ObjectDoesNotExist
from django.core.serializers.json import DjangoJSONEncoder
from django.db.models import QuerySet
from lxml import etree

from ..models import Analysis, AnalysisData, Annotation, Interaction, MetaKey, Regulation


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


def cache_result(obj, cache_name: Union[str, Path]):
    encoding = mimetypes.guess_type(cache_name)[1]
    if encoding == 'gzip':
        cache = gzip.open(cache_name, 'wb')
    else:
        cache = open(cache_name, 'wb')

    with closing(cache) as c:
        pickle.dump(obj, c, protocol=pickle.HIGHEST_PROTOCOL)


def read_cached_result(cache_name: Union[str, Path]):
    encoding = mimetypes.guess_type(cache_name)[1]
    if encoding == 'gzip':
        cache = gzip.open(cache_name, 'rb')
    else:
        cache = open(cache_name, 'rb')

    with closing(cache) as c:
        return pickle.load(c)


def convert_float(s) -> Optional[float]:
    try:
        return float(s)
    except (ValueError, TypeError):
        return None


def metadata_to_dict(df: pd.DataFrame) -> Dict[str, Any]:
    df.reset_index(inplace=True)
    meta_columns = [{"id": column, "name": column, "field": column} for column in df.columns]

    df = df.where(df.notna(), None)

    return {'columns': meta_columns,
            'data': df.to_dict(orient='index')}


nan_regex = re.compile(r'^n/?an?$', flags=re.I)

analysis_datatypes = ['ANALYSIS_METHOD', 'ANALYSIS_CUTOFF']
searchable = ['TRANSCRIPTION_FACTOR_ID', 'EXPERIMENT_TYPE', 'EXPERIMENTER', 'DATE', 'TECHNOLOGY', 'ANALYSIS_METHOD',
              'ANALYSIS_CUTOFF', 'EDGE_TYPE', 'GENOTYPE', 'DATA_SOURCE', 'TREATMENTS', 'CONTROL', 'TISSUE/SAMPLE']


def process_meta_file(f) -> pd.Series:
    metadata = pd.Series(f.readlines())

    metadata = (metadata
                .str.split(':', 1, True)
                .apply(lambda x: x.str.strip())
                .replace([r'', nan_regex], np.nan, regex=True)
                .dropna(subset=[0])
                .fillna('')
                .set_index(0, verify_integrity=True)
                .squeeze())

    metadata.index = metadata.index.str.upper().str.replace(' ', '_')

    metadata[metadata.index.str.contains(r'_?DATE$')] = pd.to_datetime(
        metadata[metadata.index.str.contains(r'_?DATE$')], infer_datetime_format=True).dt.strftime('%Y-%m-%d')

    return metadata


def process_data(f) -> Tuple[pd.DataFrame, bool]:
    data = pd.read_csv(f, header=0, na_values=['#DIV/0!', '#N/A!', '#NAME?', '#NULL!', '#NUM!', '#REF!', '#VALUE!'])
    data = data.dropna(axis=0, how='all').dropna(axis=1, how='all')

    if data.shape[1] == 3:
        data.columns = ['gene_id', 'log2fc', 'pvalue']

        data['log2fc'] = data['log2fc'].mask(np.isneginf(data['log2fc']), -sys.float_info.max)
        data['log2fc'] = data['log2fc'].mask(np.isposinf(data['log2fc']), sys.float_info.max)

        return data, True
    elif data.shape[1] == 1:
        data.columns = ['gene_id']
        return data, False
    else:
        raise ValueError(
            "Malformed Data. Must have 1 gene id column, optionally accompanied by 2 columns, log2 fold change and "
            "adjusted p-value.")


def get_exp_type(metadata: pd.Series) -> str:
    if metadata['EXPERIMENT_TYPE'] == 'EXPRESSION':
        return metadata['EXPRESSION_TYPE']
    elif metadata['EXPERIMENT_TYPE'] == 'BINDING':
        return metadata['BINDING_TYPE']
    else:
        raise ValueError('Invalid EXPERIMENT_TYPE')


def insert_data(data_file, metadata_file):
    data, has_pvals = process_data(data_file)

    with open(metadata_file) as m:
        metadata = process_meta_file(m)

    try:
        tf = Annotation.objects.get(gene_id=metadata['TRANSCRIPTION_FACTOR_ID'])
    except Annotation.DoesNotExist:
        raise ValueError('Transcription Factor ID {} does not exist.'.format(metadata['TRANSCRIPTION_FACTOR_ID']))

    # Insert Analysis
    analysis = Analysis(tf=tf)
    analysis.save()

    meta_keys = [MetaKey.objects.get_or_create(name=n, defaults={'searchable': n in searchable})
                 for n in metadata.index]

    meta_key_frame = pd.DataFrame(((m.id, m.name, m.searchable, c) for m, c in meta_keys),
                                  columns=['id', 'name', 'searchable', 'created'])
    meta_key_frame = meta_key_frame.set_index('name')

    AnalysisData.objects.bulk_create(
        [AnalysisData(analysis=analysis, key_id=meta_key_frame.at[key, 'id'], value=val)
         for key, val in metadata.iteritems()])

    anno = pd.DataFrame(Annotation.objects.filter(
        gene_id__in=data.iloc[:, 0]
    ).values_list('gene_id', 'id', named=True).iterator())

    data = data.merge(anno, on='gene_id')

    Interaction.objects.bulk_create(
        Interaction(
            analysis=analysis,
            target_id=row.id
        ) for row in data.itertuples()
    )

    if has_pvals:
        Regulation.objects.bulk_create(
            Regulation(
                analysis=analysis,
                foldchange=row.log2fc,
                p_value=row.pvalue,
                target_id=row.id
            ) for row in data.itertuples(index=False)
        )


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
        name, criterion = m.groups('')
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
