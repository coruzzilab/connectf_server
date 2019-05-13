import os

from django.apps import AppConfig

from .models import Annotation


class QuerytgdbConfig(AppConfig):
    name = 'querytgdb'

    def ready(self):
        from django.core.checks import register, Error, Warning
        from django.conf import settings

        @register()
        def check_annotations(app_configs, **kwargs):
            errors = []

            if not Annotation.objects.exists():
                errors.append(Warning('No annotations loaded.'))

            return errors

        @register()
        def check_data_files(app_configs, **kwargs):
            errors = []

            if not os.path.isfile(settings.MOTIF_ANNOTATION):
                errors.append(
                    Error('Motif annotations not found.', id='querytgdb.E001')
                )

            if not os.path.isfile(settings.MOTIF_CLUSTER):
                errors.append(
                    Error('Motif cluster info not found.', id='querytgdb.E002')
                )

            if not os.path.isdir(settings.GENE_LISTS):
                errors.append(
                    Warning('Gene list folder not found.', id='querytgdb.W001')
                )

            if not os.path.isdir(settings.TARGET_NETWORKS):
                errors.append(
                    Warning('Network folder not found.', id='querytgdb.W002')
                )

            return errors
