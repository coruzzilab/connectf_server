from django.core.management.base import BaseCommand, CommandError
from django.db.transaction import atomic

from ...utils import insert_data


class MetaDataException(CommandError):
    def __init__(self, key, *args, **kwargs):
        super().__init__("Please provide {} in metadata file".format(key), *args, **kwargs)


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('data', help='Input file for gene list')
        parser.add_argument('metadata', help='Input file for metadata of an experiment')

    def handle(self, *args, **options):
        with atomic():
            try:
                insert_data(options['data'], options['metadata'])
            except ValueError as e:
                raise CommandError(e) from e
            except KeyError as e:
                raise MetaDataException(e) from e
