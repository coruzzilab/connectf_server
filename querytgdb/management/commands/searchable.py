from django.core.management.base import BaseCommand, CommandParser
from django.db.transaction import atomic

from ...models import MetaKey


class Command(BaseCommand):
    """
    Mark metadata fields as searchable.
    """

    def add_arguments(self, parser: CommandParser):
        group = parser.add_mutually_exclusive_group()
        group.add_argument('-l', '--list', help="List metadata fields. Searchable fields marked with '*'.",
                           action="store_true")
        group.add_argument('-a', '--add', nargs='+', help='make fields searchable')
        group.add_argument('-r', '--remove', nargs='+', help='make fields not searchable')

    def handle(self, *args, **options):
        with atomic():
            if options['list']:
                self.stdout.writelines(
                    [('* ' if s else '  ') + n for n, s in MetaKey.objects.values_list('name', 'searchable')])

            elif options['add']:
                MetaKey.objects.filter(name__in=options['add']).update(searchable=True)

            elif options['remove']:
                MetaKey.objects.filter(name__in=options['remove']).update(searchable=False)
