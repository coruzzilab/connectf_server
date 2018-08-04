import glob
import os.path as path
import re

import pandas as pd
from django.core.management.base import BaseCommand, CommandError
from django.db.transaction import atomic

from ..commands.insert_tgdb import MetaDataException
from ...utils import insert_data


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('folder', help='path to tgdb data folder')
        parser.add_argument('-a', '--auto-naming', help='auto rename duplicate IDs', action='store_true')

    def handle(self, *args, **options):
        data = pd.DataFrame(glob.glob(path.join(options['folder'], 'data/*.csv')))
        data[1] = data[0].apply(path.basename).str.extract(r'^(.+)\.csv$')

        meta = pd.DataFrame(glob.glob(path.join(options['folder'], 'metadata/*.txt')))
        meta[1] = meta[0].apply(path.basename).str.extract(r'^(.+)\.txt$')

        df = data.merge(meta, on=1, how='outer')

        with atomic():
            try:
                for row in df.itertuples():
                    insert_data(row[1], row[3], auto_naming=options['auto_naming'])
            except ValueError as e:
                raise CommandError(e) from e
            except KeyError as e:
                raise MetaDataException(e) from e
