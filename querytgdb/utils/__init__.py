import base64
import gzip
import mimetypes
import os.path as path
import pickle
import random
import re
import string
from contextlib import closing
from operator import attrgetter, itemgetter
from pathlib import Path
from typing import Any, Dict, Optional, Tuple, Union
from uuid import UUID

import numpy as np
import pandas as pd
from django.conf import settings
from django.core.management.base import CommandError
from django.core.serializers.json import DjangoJSONEncoder
from lxml import etree

from ..models import Analysis, AnalysisData, Annotation, Edge, Experiment, ExperimentData, Interaction, Regulation


class PandasJSONEncoder(DjangoJSONEncoder):
    def default(self, o):
        if isinstance(o, np.number):
            if np.isinf(o):
                return float('inf')
            if np.isnan(o):
                return None
            if isinstance(o, np.integer):
                return int(o)
            if isinstance(o, np.floating):
                return float(o)
        if isinstance(o, pd.Index):
            return o.tolist()
        if isinstance(o, np.ndarray):
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


def rand_string(length):
    return ''.join(random.choices(string.ascii_letters + string.digits + '_', k=length))


nan_regex = re.compile(r'^nan?$', flags=re.I)

analysis_datatypes = ['ANALYSIS_ID', 'ANALYSIS_METHOD', 'ANALYSIS_CUTOFF', 'ANALYSIS_COMMAND', 'ANALYSIS_NOTES',
                      'ANALYSIS_BATCHEFFECT']

experiment_types = ["EXPERIMENT", "EXPERIMENT_TYPE", "EXPRESSION_TYPE", "BINDING_TYPE"]


def process_meta_file(f) -> pd.Series:
    metadata = pd.Series(f.readlines())

    metadata = (metadata
                .str.split(':', 1, True)
                .apply(lambda x: x.str.strip())
                .replace([r'', nan_regex], np.nan, regex=True)
                .set_index(0, verify_integrity=True))[1].dropna()

    metadata.index = metadata.index.str.upper().str.replace(' ', '_')

    metadata[metadata.index.str.endswith('_DATE')] = pd.to_datetime(
        metadata[metadata.index.str.endswith('_DATE')], infer_datetime_format=True).dt.strftime('%Y-%m-%d')
    metadata[metadata.index.str.endswith('_ID')] = metadata[metadata.index.str.endswith('_ID')].str.upper()
    metadata[metadata.index.isin(experiment_types)] = metadata[metadata.index.isin(experiment_types)].str.upper()

    return metadata


def process_data(f) -> Tuple[pd.DataFrame, bool]:
    data = pd.read_csv(f, header=0, index_col=0)

    has_pvals = True

    try:
        # with fc pvalue induce/repress edge
        data.iloc[:, [2, 3]] = data.iloc[:, [2, 3]].apply(lambda x: x.str.upper())
    except IndexError:
        # edge only
        data.iloc[:, 0] = data.iloc[:, 0].str.upper()
        has_pvals = False

    return data, has_pvals


def get_exp_type(metadata: pd.Series) -> str:
    if metadata['EXPERIMENT_TYPE'] == 'EXPRESSION':
        return metadata['EXPRESSION_TYPE']
    elif metadata['EXPERIMENT_TYPE'] == 'BINDING':
        return metadata['BINDING_TYPE']
    else:
        raise ValueError('Invalid EXPERIMENT_TYPE')


def insert_data(data_file, metadata_file, auto_naming=False):
    data, has_pvals = process_data(data_file)

    with open(metadata_file) as m:
        metadata = process_meta_file(m)

    try:
        tf = Annotation.objects.get(gene_id=metadata['TRANSCRIPTION_FACTOR_ID'])
    except Annotation.DoesNotExist:
        raise ValueError('Transcription Factor ID does not exist.')

    try:
        edge_prefix = metadata['EXPERIMENT'] + ':' + get_exp_type(metadata) + ':'
    except ValueError as e:
        raise CommandError(e) from e

    edge_prefix = edge_prefix.upper()

    try:
        exp_id = metadata['EXPERIMENT_ID']
    except KeyError:
        tf_name, experimenter, date_ = itemgetter('TRANSCRIPTION_FACTOR_ID', 'EXPERIMENTER',
                                                  'EXPERIMENT_DATE')(metadata)

        exp_id = tf_name + '_' + ''.join(
            map(itemgetter(0), re.split(r'\s+|_', experimenter))) + '{:%m%d%y}_'.format(pd.to_datetime(date_))

        exp_id += get_exp_type(metadata)

    exp_id = exp_id.upper()

    # Insert Experiment
    experiment, exp_created = Experiment.objects.get_or_create(name=exp_id, tf=tf)

    try:
        analysis_id = metadata['ANALYSIS_ID']
    except KeyError:
        analysis_id = metadata['TRANSCRIPTION_FACTOR_ID'] + '_' + metadata[
            'ANALYSIS_METHOD'] + '_' + rand_string(10)

    if Analysis.objects.filter(experiment=experiment, name=analysis_id).exists() and auto_naming:
        analysis_id += '_' + rand_string(10)

    analysis = Analysis(experiment=experiment, name=analysis_id)
    analysis.save()

    if exp_created:
        # Insert experiment metadata if new experiment
        ExperimentData.objects.bulk_create(
            [ExperimentData(experiment=experiment, key=key, value=val) for key, val in
             metadata[~metadata.index.isin(analysis_datatypes)].iteritems()])

    AnalysisData.objects.bulk_create(
        [AnalysisData(analysis=analysis, key=key, value=val) for key, val in
         metadata[metadata.index.isin(analysis_datatypes)].iteritems()])

    if has_pvals:
        edges = (edge_prefix +
                 data.iloc[:, 2].str.cat(data.iloc[:, 3], ':'))
    else:
        edges = edge_prefix + data.iloc[:, 0]

    edges.name = "edge"
    edges.index.name = None

    data = pd.concat([data, edges], axis=1)

    edge_objs = map(itemgetter(0), (Edge.objects.get_or_create(name=e) for e in edges.unique()))
    edge_frame = pd.DataFrame.from_records(
        map(attrgetter('id', 'name'), edge_objs), columns=['edge_id', 'edge'])
    data = data.reset_index().merge(edge_frame, on='edge', how='left')

    anno = pd.DataFrame(Annotation.objects.filter(
        gene_id__in=data['index']
    ).values_list('gene_id', 'id', named=True).iterator())

    data = data.merge(anno, left_on='index', right_on='gene_id')

    Interaction.objects.bulk_create(
        Interaction(
            analysis=analysis,
            target_id=row.id,
            edge_id=row.edge_id
        ) for row in data.itertuples()
    )

    if has_pvals:
        Regulation.objects.bulk_create(
            Regulation(
                analysis=analysis,
                foldchange=row[1],
                p_value=row[2],
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
