from argparse import ArgumentParser
from operator import attrgetter, itemgetter

import pandas as pd
from django.core.management.base import BaseCommand
from django.db.transaction import atomic

from querytgdb.models import Annotation, EdgeData, EdgeType
from ...utils.sif import get_network


class Command(BaseCommand):
    help = "Adds edge properties to annotate gene interactions in the database, DAP, DAP_amp, etc."

    def add_arguments(self, parser: ArgumentParser):
        parser.add_argument('file', help='edge property file', nargs='?')
        parser.add_argument('--clear', help='clear current edges', action='store_true')

    def handle(self, *args, **options):
        with atomic():
            if options['clear']:
                EdgeData.objects.all().delete()
                EdgeType.objects.all().delete()

            if 'file' in options:
                with open(options["file"]) as f:
                    g = get_network(f)
                df = pd.DataFrame(iter(g.edges(keys=True)))
                df.columns = ['source', 'target', 'edge']
                df = df.drop_duplicates()

                edges = pd.DataFrame.from_records(map(attrgetter('id', 'name'),
                                                      map(itemgetter(0),
                                                          (EdgeType.objects.get_or_create(name=e) for e in
                                                           df['edge'].unique()))),
                                                  columns=['edge_id', 'edge'])

                anno = pd.DataFrame(Annotation.objects.values_list('id', 'gene_id', named=True).iterator())

                df = (df
                      .merge(edges, on='edge')
                      .merge(anno, left_on='source', right_on='gene_id')
                      .merge(anno, left_on='target', right_on='gene_id'))

                EdgeData.objects.bulk_create(
                    (EdgeData(
                        type_id=e,
                        tf_id=s,
                        target_id=t
                    ) for e, s, t in df[['edge_id', 'id_x', 'id_y']].itertuples(index=False, name=None)),
                    batch_size=1000
                )
