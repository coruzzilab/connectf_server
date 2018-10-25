from django.core.management.base import BaseCommand, CommandParser
from django.db.transaction import atomic

from ...models import Analysis


class Command(BaseCommand):
    help = "removing experiments from the database"

    def add_arguments(self, parser: CommandParser):
        parser.add_argument("-d", "--dry-run", help="print experiments to be deleted", action="store_true")
        parser.add_argument("--all", help="Delete all data", action="store_true")

        subparsers = parser.add_subparsers(dest="command")

        parser_analysis = subparsers.add_parser('analysis', help="remove data by analysis ID")
        parser_analysis.add_argument("analysis_id", nargs='+', help="analysis IDs to remove", type=int)

        parser_meta = subparsers.add_parser('metadata', help="remove data by metadata fields")
        parser_meta.add_argument("key", help="metadata key")
        parser_meta.add_argument("value", help="metadata value")

    def handle(self, *args, **options):
        dry_run = options["dry_run"]
        command = options["command"]

        with atomic():
            if command == "analysis":
                queryset = Analysis.objects.filter(pk__in=options["analysis_id"])

                for analysis in queryset:
                    self.stdout.write("Analysis: {}".format(analysis.name), ending="\n")

                if not dry_run:
                    self.stdout.write("deleting analysis {0[analysis_id]}".format(options), ending="\n")
                    queryset.delete()

            elif command == "metadata":
                analyses = Analysis.objects.filter(
                    analysisdata__key__name__iexact=options["key"],
                    analysisdata__value__icontains=options["value"]
                )

                for analysis in analyses:
                    self.stdout.write("Analysis: {}".format(analysis.name),
                                      ending="\n")

                if not dry_run:
                    self.stdout.write("deleting analyses with key: {0[key]} value: {0[value]}".format(options),
                                      ending="\n")
                    analyses.delete()

            elif command is None and options["all"]:
                analyses = Analysis.objects.all()

                self.stdout.write("deleting everything...", ending="\n")
                analyses.delete()

            else:
                self.stdout.write("Nothing removed.", ending="\n")
