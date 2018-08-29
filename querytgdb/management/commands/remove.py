from django.core.management.base import BaseCommand, CommandParser
from django.db.transaction import atomic

from ...models import Analysis, Experiment


class Command(BaseCommand):
    help = "removing experiments from the database"

    def add_arguments(self, parser: CommandParser):
        parser.add_argument("-d", "--dry-run", help="print experiments to be deleted", action="store_true")

        subparsers = parser.add_subparsers(dest="command")

        parser_exp = subparsers.add_parser('experiment', help="remove data by experiment ID")
        parser_exp.add_argument("exp_id", nargs='+', help="experiment IDs to remove")

        parser_analysis = subparsers.add_parser('analysis', help="remove data by analysis ID")
        parser_analysis.add_argument("analysis_id", nargs='+', help="analysis IDs to remove")

        parser_meta = subparsers.add_parser('metadata', help="remove data by metadata fields")
        parser_meta.add_argument("key", help="metadata key")
        parser_meta.add_argument("value", help="metadata value")

    def handle(self, *args, **options):
        dry_run = options["dry_run"]
        command = options["command"]

        with atomic():
            if not dry_run:
                self.stdout.write("deleting...", ending="\n")

            if command == "experiment":
                queryset = Experiment.objects.filter(name__in=options["exp_id"]).prefetch_related("analysis_set")

                for exp in queryset:
                    for analysis in exp.analysis_set.all():
                        self.stdout.write("Experiment: {} Analysis: {}".format(exp.name, analysis.name), ending="\n")

                if not dry_run:
                    queryset.delete()

            elif command == "analysis":
                queryset = Analysis.objects.filter(name__in=options["analysis_id"]).prefetch_related("experiment")

                for analysis in queryset:
                    self.stdout.write("Experiment: {} Analysis: {}".format(analysis.experiment.name, analysis.name),
                                      ending="\n")

                if not dry_run:
                    queryset.delete()

            elif command == "metadata":
                analyses = Analysis.objects.filter(
                    analysisdata__key__iexact=options["key"],
                    analysisdata__value__icontains=options["value"]
                ).prefetch_related("experiment")

                experiments = Experiment.objects.filter(
                    experimentdata__key__iexact=options["key"],
                    experimentdata__value__icontains=options["value"]
                ).prefetch_related("analysis_set")

                for analysis in analyses:
                    self.stdout.write("Experiment: {} Analysis: {}".format(analysis.experiment.name, analysis.name),
                                      ending="\n")

                for exp in experiments:
                    for analysis in exp.analysis_set.all():
                        self.stdout.write("Experiment: {} Analysis: {}".format(exp.name, analysis.name), ending="\n")

                if not dry_run:
                    analyses.delete()
                    experiments.delete()
