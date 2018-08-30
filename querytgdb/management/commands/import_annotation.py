from argparse import ArgumentParser
from operator import itemgetter

import pandas as pd
from django.core.management.base import BaseCommand
from django.db.transaction import atomic

from ...models import Annotation


class Command(BaseCommand):
    help = "Import or update Annotations used to annotate results"

    def add_arguments(self, parser: ArgumentParser):
        group = parser.add_mutually_exclusive_group(required=True)

        group.add_argument("-i", "--input", help="annotation file")
        group.add_argument("-o", "--output", help="export current annotations")

    def handle(self, *args, **options):
        infile, outfile = itemgetter("input", "output")(options)

        anno = pd.DataFrame(Annotation.objects.values_list(named=True).iterator())

        if outfile:
            anno[["gene_id", "name", "fullname", "gene_type", "gene_family"]].to_csv(outfile, index=False)

        elif infile:
            anno = anno.set_index('gene_id').fillna('')

            in_anno = pd.read_csv(infile)
            in_anno.columns = ["gene_id", "name", "fullname", "gene_type", "gene_family"]
            in_anno = in_anno.set_index('gene_id').fillna('')

            changed = (in_anno.loc[anno.index, ["name", "fullname", "gene_type", "gene_family"]] != anno[
                ["name", "fullname", "gene_type", "gene_family"]]).any(axis=1)

            to_update = pd.concat([
                anno['id'],
                in_anno.loc[changed[changed].index, :]
            ], axis=1, join='inner').reset_index()

            new_anno = in_anno.loc[~in_anno.index.isin(anno.index), :].reset_index()

            with atomic():
                for a in (Annotation(**row._asdict()) for row in to_update.itertuples(index=False)):
                    a.save()

                Annotation.objects.bulk_create(
                    (Annotation(**row._asdict()) for row in new_anno.itertuples(index=False)))
