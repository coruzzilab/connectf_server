import re
from operator import attrgetter, itemgetter
from typing import Tuple

import numpy as np
import pandas as pd
from django.core.management.base import BaseCommand, CommandError
from django.db.models import Max
from django.db.transaction import atomic

from ...models import Analysis, AnalysisIddata, Annotation, Edges, Interactions, MetaIddata, Metadata, ReferenceId, \
    Regulation, TargetDBTF
from ...utils import rand_string

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


class MetaDataException(CommandError):
    def __init__(self, key, *args, **kwargs):
        super().__init__("Please provide {} in metadata file".format(key), *args, **kwargs)


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('data', help='Input file for gene list')
        parser.add_argument('metadata', help='Input file for metadata of an experiment')

    def handle(self, *args, **options):
        with atomic():
            data, has_pvals = process_data(options['data'])

            try:
                with open(options['metadata']) as m:
                    metadata = process_meta_file(m)
            except ValueError as e:
                raise CommandError(e) from e

            try:
                tf, created = TargetDBTF.objects.get_or_create(db_tf_agi=metadata['TRANSCRIPTION_FACTOR_ID'])
            except KeyError as e:
                raise MetaDataException("TRANSCRIPTION_FACTOR_ID") from e

            try:
                try:
                    edge_prefix = metadata['EXPERIMENT'] + ':' + get_exp_type(metadata) + ':'
                except ValueError as e:
                    raise CommandError(e) from e

                edge_prefix = edge_prefix.upper()

            except KeyError as e:
                raise MetaDataException(e) from e

            try:
                exp_id = metadata['EXPERIMENT_ID']
            except KeyError:
                tf_name, experimenter, date_ = itemgetter('TRANSCRIPTION_FACTOR_ID', 'EXPERIMENTER',
                                                          'EXPERIMENT_DATE')(metadata)

                exp_id = tf_name + '_' + ''.join(
                    map(itemgetter(0), re.split(r'\s+|_', experimenter))) + '{:%m%d%y}_'.format(pd.to_datetime(date_))

                try:
                    exp_id += get_exp_type(metadata)
                except ValueError as e:
                    raise CommandError(e) from e

            exp_id = exp_id.upper()

            # Insert Metadata
            try:
                meta = Metadata.objects.get(meta_fullid=exp_id)
            except Metadata.DoesNotExist:
                new_meta_id = Metadata.objects.aggregate(new_meta_id=Max('meta_id') + 1)['new_meta_id']
                if new_meta_id is None:
                    new_meta_id = 1

                meta = Metadata(meta_id=new_meta_id, meta_fullid=exp_id)
                meta.save()
            except KeyError as e:
                raise MetaDataException(e) from e

            # Insert experiment metadata
            try:
                try:
                    analysis_id = metadata['ANALYSIS_ID']
                except KeyError:
                    analysis_id = metadata['TRANSCRIPTION_FACTOR_ID'] + '_' + metadata[
                        'ANALYSIS_METHOD'] + '_' + rand_string(10)
            except KeyError as e:
                raise MetaDataException(e) from e

            try:
                if ReferenceId.objects.filter(meta_id=meta,
                                              analysis_id__analysis_fullid=analysis_id).exists():
                    # dev purposes only
                    # ReferenceId.objects.filter(meta_id=meta,
                    #                            analysis_id__analysis_fullid=analysis_id).delete()
                    # Analysis.objects.filter(analysis_fullid=analysis_id).delete()

                    raise CommandError('ANALYSIS_ID already exists')
            except KeyError as e:
                raise MetaDataException('ANALYSIS_ID') from e

            new_analysis_id = Analysis.objects.aggregate(new_analysis_id=Max('analysis_id') + 1)['new_analysis_id']
            if new_analysis_id is None:
                new_analysis_id = 1

            analysis = Analysis(analysis_fullid=analysis_id, analysis_id=new_analysis_id)
            analysis.save()

            ref = ReferenceId(meta_id=meta, analysis_id=analysis)
            ref.save()

            MetaIddata.objects.bulk_create(
                [MetaIddata(meta_id=meta, meta_type=key, meta_value=val) for key, val in
                 metadata[~metadata.index.isin(analysis_datatypes)].iteritems()])

            AnalysisIddata.objects.bulk_create(
                [AnalysisIddata(analysis_id=analysis, analysis_type=key, analysis_value=val) for key, val in
                 metadata[metadata.index.isin(analysis_datatypes)].iteritems()])

            try:
                if has_pvals:
                    edges = (edge_prefix +
                             data.iloc[:, 2].str.cat(data.iloc[:, 3], ':'))
                else:
                    edges = edge_prefix + data.iloc[:, 0]

                edges.name = "edge"
                edges.index.name = None

                data = pd.concat([data, edges], axis=1)

                edge_objs = map(itemgetter(0), (Edges.objects.get_or_create(edge_name=e) for e in edges.unique()))
                edge_frame = pd.DataFrame.from_records(
                    map(attrgetter('edge_id', 'edge_name'), edge_objs), columns=['edge_id', 'edge'])
                data = data.reset_index().merge(edge_frame, on='edge', how='left')

                anno = pd.DataFrame(Annotation.objects.filter(
                    agi_id__in=data['index']
                ).values_list('agi_id', 'ath_id', named=True).iterator())

                data = data.merge(anno, left_on='index', right_on='agi_id')

                Interactions.objects.bulk_create(
                    Interactions(
                        ref_id=ref,
                        db_tf_id=tf,
                        target_id_id=row.ath_id,
                        edge_id_id=row.edge_id
                    ) for row in data.itertuples()
                )

                if has_pvals:
                    Regulation.objects.bulk_create(
                        Regulation(
                            ref_id=ref,
                            ath_id_id=row.ath_id,
                            foldchange=row[1],
                            pvalue=row[2]
                        ) for row in data.astype(str).itertuples(index=False)
                    )

            except KeyError as e:
                raise MetaDataException(e) from e
