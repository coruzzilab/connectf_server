import os

from django.core.management.base import BaseCommand, CommandParser
from django.db.transaction import atomic

from querytgdb.models import EdgeData, EdgeType
from querytgdb.utils.insert_data import import_additional_edges


class Command(BaseCommand):
    help = "Adds edge properties to annotate gene interactions in the database, DAP, DAP_amp, etc."

    def add_arguments(self, parser: CommandParser):
        group = parser.add_mutually_exclusive_group()

        group2 = group.add_argument_group(title='file import')
        group2.add_argument('file', help='edge property file', nargs='?')
        group2.add_argument('-f', '--format', help='file format', type=str)
        group2.add_argument('-U', '--undirected', help='treat edges as undirected', action='store_true')

        group.add_argument('--clear-all', help='clear all edges', action='store_true')
        group.add_argument('--clear', help='clear edge', type=str)
        group.add_argument('-l', '--list', action='store_true')

    def handle(self, *args, **options):
        with atomic():
            if options['clear_all']:
                EdgeData.objects.all().delete()
                EdgeType.objects.all().delete()

            if options['list']:
                for edge in EdgeType.objects.all():
                    print(edge.name)

            if options['clear'] is not None:
                qs = EdgeType.objects.filter(name=options['clear'])
                if qs.exists():
                    print(f"Deleting: {options['clear']}")
                    qs.delete()

            if options['file'] is not None:
                name, ext = os.path.splitext(options['file'])

                import_additional_edges(options["file"], sif=(ext == '.sif'), directional=(not options['undirected']))
