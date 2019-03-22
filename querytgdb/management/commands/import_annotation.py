from operator import itemgetter

import pandas as pd
from django.core.management.base import BaseCommand, CommandParser

from querytgdb.utils.insert_data import import_annotations
from ...models import Annotation


class Command(BaseCommand):
    help = "Import or update Annotations used to annotate results"

    def add_arguments(self, parser: CommandParser):
        parser.add_argument("--delete", help="delete data from database if it doesn't appear in import",
                            action="store_true")
        parser.add_argument("-d", "--dry-run", help="print changes only", action="store_true")
        group = parser.add_mutually_exclusive_group(required=True)

        group.add_argument("-i", "--input", help="annotation file")
        group.add_argument("-o", "--output", help="export current annotations")

    def handle(self, *args, **options):
        infile, outfile = itemgetter("input", "output")(options)

        if outfile:
            anno = pd.DataFrame(Annotation.objects.values_list(named=True).iterator())
            anno[["gene_id", "name", "fullname", "gene_type", "gene_family"]].to_csv(outfile, index=False)

        elif infile:
            import_annotations(infile, dry_run=options['dry_run'], delete_existing=options['delete'])
