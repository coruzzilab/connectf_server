import glob
import os.path as path

import pandas as pd
from django.core.management.base import BaseCommand, CommandError
from django.db.transaction import atomic

from querytgdb.utils.insert_data import insert_data


def get_file_name(x):
    return path.splitext(path.basename(x))[0]


class MetaDataException(CommandError):
    def __init__(self, key, *args, **kwargs):
        super().__init__("Please provide {} in metadata file".format(key), *args, **kwargs)


class Command(BaseCommand):
    help = "Loads data into the application."

    def add_arguments(self, parser):
        parser.add_argument('data', help='Data file or folder for an experiment')
        parser.add_argument('metadata', help='Metadata file or folder for an experiment')
        parser.add_argument('--sep', type=str, help="specify separator for csv. [default: ',']", default=',')
        parser.add_argument('-d', '--dry-run', action='store_true', help="dry run")

    def handle(self, *args, **options):
        with atomic():
            try:
                if path.isdir(options["data"]):
                    if not path.isdir(options["metadata"]):
                        raise CommandError("Both data and metadata should be folders.")

                    data = pd.DataFrame(glob.glob(path.join(options['data'], '*')))
                    data[1] = data[0].apply(get_file_name)

                    meta = pd.DataFrame(glob.glob(path.join(options['metadata'], '*')))
                    meta[1] = meta[0].apply(get_file_name)

                    df = data.merge(meta, on=1, how='inner', validate='one_to_one')

                    for row in df.itertuples():
                        insert_data(row[1], row[3], sep=options["sep"], dry_run=options['dry_run'])
                else:
                    insert_data(options['data'], options['metadata'], sep=options["sep"], dry_run=options['dry_run'])
            except ValueError as e:
                raise CommandError(e) from e
            except KeyError as e:
                raise MetaDataException(e) from e
