import re
import sys
from typing import Tuple

import numpy as np
import pandas as pd

from querytgdb.models import Analysis, AnalysisData, Annotation, Interaction, MetaKey, Regulation

nan_regex = re.compile(r'^n/?an?$', flags=re.I)
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


def insert_data(data_file, metadata_file, sep=','):
    data, has_pvals = process_data(data_file, sep=sep)

    with open(metadata_file) as m:
        metadata = process_meta_file(m)

    try:
        tf = Annotation.objects.get(gene_id=metadata['TRANSCRIPTION_FACTOR_ID'])
    except Annotation.DoesNotExist:
        raise ValueError('Transcription Factor ID {} does not exist.'.format(metadata['TRANSCRIPTION_FACTOR_ID']))

    if not metadata.index.contains('EDGE_TYPE'):
        raise ValueError('Please assign an EDGE_TYPE to the metadata.')

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
