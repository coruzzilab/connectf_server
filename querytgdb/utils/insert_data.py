import logging
import re
import sys
import warnings
from operator import attrgetter, itemgetter
from typing import TextIO, Tuple

import numpy as np
import pandas as pd
from django.db.transaction import atomic
from django.db.utils import IntegrityError

from querytgdb.models import Analysis, AnalysisData, Annotation, EdgeData, EdgeType, Interaction, MetaKey, Regulation, ImportHistory
from querytgdb.utils.sif import get_network

logger = logging.getLogger(__name__)

nan_regex = re.compile(r'^n/?an?$', flags=re.I)


class UnkownGeneWarning(Warning):
    pass


def process_meta_file(f: TextIO) -> pd.DataFrame:
    metadata = pd.Series(f.readlines())

    metadata = (metadata
                .str.split(':', n=1, expand=True)
                .apply(lambda x: x.str.strip())
                .replace([r'', nan_regex], np.nan, regex=True)
                .dropna(subset=[0])
                .fillna('')
                .set_index(0, verify_integrity=True))

    metadata.index = metadata.index.str.upper().str.replace(' ', '_')
    metadata['special'] = metadata.index.str.extract(r'^([^a-z])', flags=re.I, expand=False).fillna('')
    metadata.index = metadata.index.str.replace(r'^([^a-z])', '', flags=re.I, regex=True)
    metadata.columns = ['data', 'special']

    date_rows = metadata.index.str.contains(r'_?DATE$')
    try:
        metadata.loc[date_rows, 'data'] = pd.to_datetime(
            metadata.loc[date_rows, 'data'], infer_datetime_format=True).dt.strftime('%Y-%m-%d')
    except ValueError as e:
        pass

    return metadata


def process_data(f, sep=',') -> Tuple[pd.DataFrame, bool]:
    data = pd.read_csv(f, header=0, sep=sep,
                       na_values=['#DIV/0!', '#N/A!', '#NAME?', '#NULL!', '#NUM!', '#REF!', '#VALUE!'])
    data = data.dropna(axis=0, how='all').dropna(axis=1, how='all')

    cols = data.shape[1]

    if cols == 3:
        data.columns = ['gene_id', 'log2fc', 'pvalue']

        data['log2fc'] = data['log2fc'].mask(np.isneginf(data['log2fc']), -sys.float_info.max)
        data['log2fc'] = data['log2fc'].mask(np.isposinf(data['log2fc']), sys.float_info.max)

        return data, True
    elif cols == 1:
        data.columns = ['gene_id']
        return data, False
    else:
        raise ValueError(
            "Malformed Data. Must have 1 gene id column, optionally accompanied by 2 columns, log2 fold change and "
            "adjusted p-value.")


def insert_data(data_file, metadata_file, sep=',', dry_run=False):
    # if data_file or metadata_file are already imported, skip
    if (
        ImportHistory.objects.filter(fileName=data_file, type="data").exists() or 
        ImportHistory.objects.filter(fileName=metadata_file, type="metadata").exists()
    ):
        print(f"Skipping {data_file.split('/')[-1]} and {metadata_file.split('/')[-1]} as they are already imported.")
        return

    print("inserting {0} {1}\n".format(data_file, metadata_file))

    data, has_pvals = process_data(data_file, sep=sep)

    try:
        with open(metadata_file) as m:
            metadata = process_meta_file(m)
    except TypeError:
        metadata = process_meta_file(metadata_file)

    if not Annotation.objects.exists():
        raise ValueError("Please use 'import_annotation' to add gene annotations to the database first.")

    try:
        tf = Annotation.objects.get(gene_id=metadata.loc['TRANSCRIPTION_FACTOR_ID', 'data'])
        metadata = metadata.drop('TRANSCRIPTION_FACTOR_ID')
    except Annotation.DoesNotExist:
        raise ValueError(
            'Transcription Factor ID {} does not exist.'.format(metadata.loc['TRANSCRIPTION_FACTOR_ID', 'data']))

    if 'EDGE_TYPE' not in metadata.index:
        raise ValueError('Please assign an EDGE_TYPE to the metadata.')

    try:
        exp_type = metadata.loc['EXPERIMENT_TYPE', 'data'].upper()

        if exp_type not in ('EXPRESSION', 'BINDING'):
            raise ValueError('EXPERIMENT_TYPE should be "Expression" or "Binding".')
    except KeyError:
        raise ValueError('Please assign an EXPERIMENT_TYPE to the metadata. Should be "Expression" or "Binding".')

    if dry_run:
        return

    # Insert Analysis
    analysis = Analysis(tf=tf)
    analysis.save()

    meta_keys = [MetaKey.objects.get_or_create(name=n, defaults={'searchable': s})
                 for n, s in zip(metadata.index,
                                 metadata['special'].str.contains('*', regex=False))]

    meta_key_frame = pd.DataFrame(((m.id, m.name, m.searchable, c) for m, c in meta_keys),
                                  columns=['id', 'name', 'searchable', 'created'])
    meta_key_frame = meta_key_frame.set_index('name')

    AnalysisData.objects.bulk_create(
        [AnalysisData(analysis=analysis, key_id=meta_key_frame.at[key, 'id'], value=val)
         for key, val in metadata['data'].items()])

    anno = pd.DataFrame(Annotation.objects.filter(
        gene_id__in=data.iloc[:, 0]
    ).values_list('gene_id', 'id', named=True).iterator())

    # update import history
    ImportHistory.objects.create(fileName=data_file, type="data")
    ImportHistory.objects.create(fileName=metadata_file, type="metadata")

    # check for invalid ids
    data = data.merge(anno, on='gene_id', how='left')

    unknown_genes = data['id'].isnull()

    if unknown_genes.any():
        warnings.warn("Unkown gene ids: {}".format(", ".join(data.loc[unknown_genes, 'gene_id'].unique())),
                      UnkownGeneWarning)
        data = data[~unknown_genes]

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


def read_annotation_file(annotation_file: str) -> pd.DataFrame:
    in_anno = pd.read_csv(annotation_file, comment='#').fillna('')
    in_anno.columns = ["gene_id", "name", "fullname", "gene_type", "gene_family"]

    return in_anno


def import_annotations(annotation_file: str, dry_run: bool = False, delete_existing: bool = True):
    logger.setLevel(logging.INFO)

    anno = pd.DataFrame(Annotation.objects.values_list(named=True).iterator())
    if anno.empty:
        anno = pd.DataFrame(columns=["id", "gene_id", "name", "fullname", "gene_type", "gene_family"])

    anno = anno.set_index('gene_id').fillna('')

    in_anno = read_annotation_file(annotation_file)
    in_anno = in_anno.set_index('gene_id')

    changed = (in_anno.reindex(index=anno.index).loc[:, ["name", "fullname", "gene_type", "gene_family"]] != anno[
        ["name", "fullname", "gene_type", "gene_family"]]).any(axis=1)

    to_update = pd.concat([
        anno['id'],
        in_anno.loc[in_anno.index.isin(changed[changed].index), :]
    ], axis=1, join='inner').reset_index()

    new_anno = in_anno.loc[~in_anno.index.isin(anno.index), :].reset_index()

    to_delete = anno.loc[anno.index.difference(in_anno.index), :]

    if dry_run:
        logger.info("Update:")
        logger.info(to_update)
        logger.info("Create:")
        logger.info(new_anno)

        if delete_existing:
            logger.info("Delete:")
            logger.info(to_delete)
    else:
        with atomic():
            for a in (Annotation(**row._asdict()) for row in to_update.itertuples(index=False)):
                a.save()

            Annotation.objects.bulk_create(
                (Annotation(**row._asdict()) for row in new_anno.itertuples(index=False)))

            if delete_existing:
                Annotation.objects.filter(pk__in=to_delete['id']).delete()


def import_additional_edges(edge_file: str, sif: bool = False, directional: bool = True):
    print(f"Inserting edge file: {edge_file}\n")

    if sif:
        try:
            with open(edge_file) as f:
                g = get_network(f)
        except TypeError:
            g = get_network(edge_file)

        df = pd.DataFrame(iter(g.edges(keys=True)))
    else:
        df = pd.read_csv(edge_file)
        df = df.dropna(axis=0, how='all').dropna(axis=1, how='all')

    df.columns = ['source', 'target', 'edge']
    df = df.drop_duplicates()

    edges = pd.DataFrame.from_records(map(attrgetter('id', 'name'),
                                          map(itemgetter(0),
                                              (EdgeType.objects.get_or_create(
                                                  name=e,
                                                  directional=directional
                                              ) for e in
                                                  df['edge'].unique()))),
                                      columns=['edge_id', 'edge'])

    anno = pd.DataFrame(Annotation.objects.values_list('id', 'gene_id', named=True).iterator())

    df = (df
          .merge(edges, on='edge')
          .merge(anno, left_on='source', right_on='gene_id', how='left')
          .merge(anno, left_on='target', right_on='gene_id', how='left'))

    # warn of unknown genes
    unknown_source = df['id_x'].isnull()
    unknown_target = df['id_y'].isnull()

    unknown_edges = unknown_source | unknown_target

    if unknown_edges.any():
        warnings.warn(
            "Unknown genes: {}".format(
                ", ".join(pd.concat([df.loc[unknown_source, 'source'], df.loc[unknown_target, 'target']]).unique())
            ),
            UnkownGeneWarning
        )
        df = df[~unknown_edges]

    df = df[['edge_id', 'id_x', 'id_y']]

    if not directional:
        und_df = df.copy()
        und_df[['id_x', 'id_y']] = und_df[['id_y', 'id_x']]
        df = pd.concat([df, und_df])
        df = df.drop_duplicates()

    try:
        EdgeData.objects.bulk_create(
            (EdgeData(
                type_id=e,
                tf_id=s,
                target_id=t
            ) for e, s, t in df[['edge_id', 'id_x', 'id_y']].itertuples(index=False, name=None)),
            batch_size=1000
        )
    except IntegrityError as e:
        print(f"Duplicate entry error, skipping edge file: {edge_file.split('/')[-1]}")
