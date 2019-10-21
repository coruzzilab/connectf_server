import gzip
import mimetypes
import os.path
import shutil

from django.conf import settings
from django.core.management.base import BaseCommand, CommandParser


def gz_copy(src, dst, force=False):
    mtype, enc = mimetypes.guess_type(src)

    if not force and os.path.isfile(dst):
        raise ValueError("Destination file exists, set 'force' to True to overwrite.")

    if mtype == 'text/csv':
        if enc is None:
            with open(src, 'rb') as f, gzip.open(dst, 'wb') as g:
                shutil.copyfileobj(f, g)

        elif enc == 'gzip':
            shutil.copyfile(src, dst)
    else:
        raise ValueError("Should be .csv file")


class Command(BaseCommand):
    def add_arguments(self, parser: CommandParser):
        parser.add_argument("desc", help="motif description (.csv)", type=str)
        parser.add_argument("motif", help="gene to motif matches (.csv)", type=str)

    def handle(self, *args, **options):
        gz_copy(options['desc'], os.path.join(settings.BASE_DIR, 'data/cluster_info.csv.gz'))
        gz_copy(options['motif'], os.path.join(settings.BASE_DIR, 'data/annotated_motifs.csv.gz'))
