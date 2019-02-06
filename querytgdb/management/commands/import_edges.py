import os

from django.core.management.base import BaseCommand, CommandParser
from django.db.transaction import atomic

from querytgdb.models import EdgeData, EdgeType
from querytgdb.utils.insert_data import import_additional_edges


class Command(BaseCommand):
    help = "Adds edge properties to annotate gene interactions in the database, DAP, DAP_amp, etc."

    def add_arguments(self, parser: CommandParser):
        parser.add_argument('file', help='edge property file', nargs='?')
        parser.add_argument('-f', '--format', help='file format', type=str)
        parser.add_argument('-U', '--undirected', help='treat edges as undirected', action='store_true')
        parser.add_argument('--clear', help='clear current edges', action='store_true')

    def handle(self, *args, **options):
        with atomic():
            if options['clear']:
                EdgeData.objects.all().delete()
                EdgeType.objects.all().delete()

            if 'file' in options:
                name, ext = os.path.splitext(options['file'])

                import_additional_edges(options["file"], sif=(ext == '.sif'), directional=(not options['undirected']))
