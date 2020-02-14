import csv
import tempfile

import pandas as pd
from django.core.management.base import BaseCommand, CommandParser
from django.db import connection
from django.db.transaction import atomic

from querytgdb.models import MotifName, MotifRegion


class Command(BaseCommand):
    def add_arguments(self, parser: CommandParser):
        parser.add_argument('motifs')
        parser.add_argument('output')

    def handle(self, *args, **options):
        df = pd.read_csv(options['motifs'], header=None, na_filter=False)
        with atomic():
            MotifRegion.objects.bulk_create([MotifRegion(name=name) for name in df[1].unique()])
            MotifName.objects.bulk_create([MotifName(name=name) for name in df[2].unique()])

            regions = dict(MotifRegion.objects.values_list('name', 'id'))
            motifs = dict(MotifName.objects.values_list('name', 'id'))

            it = df.itertuples(index=False, name=None)

            with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', newline='') as f:
                writer = csv.writer(f, dialect='unix', quoting=csv.QUOTE_MINIMAL)
                for c in it:
                    writer.writerow([c[0], regions[c[1]], c[3], motifs[c[2]], regions[c[1]]])

                with connection.cursor() as cursor:
                    cursor.execute(
                        "LOAD DATA LOCAL INFILE '%s' INTO TABLE querytgdb_motifmatch FIELDS TERMINATED BY ',' "
                        "ENCLOSED BY '\"' LINES TERMINATED BY '\\n' (gene_id, count, motif_id, region_id)",
                        [f.name])

            # for chunk in iter(partial(get_chunks, it), []):
            #     MotifMatch.objects.bulk_create(
            #         [MotifMatch(
            #             gene_id=c[0],
            #             region_id=regions[c[1]],
            #             motif_id=motifs[c[2]],
            #             count=c[3]
            #         ) for c in chunk])
