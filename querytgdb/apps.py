import os

from django.apps import AppConfig


class QuerytgdbConfig(AppConfig):
    name = 'querytgdb'

    def ready(self):
        from django.core.checks import register, Error, Warning
        from django.conf import settings
        from .models import Annotation
        from django.db.utils import DatabaseError

        @register()
        def check_annotations(app_configs, **kwargs):
            errors = []

            try:
                if not Annotation.objects.exists():
                    errors.append(Warning('No annotations loaded.', id='querytgdb.W003'))
            except DatabaseError:
                pass

            return errors

        @register()
        def check_data_files(app_configs, **kwargs):
            errors = []

            if not os.path.isfile(settings.MOTIF_ANNOTATION):
                errors.append(
                    Error('Motif annotations not found.', id='querytgdb.E001')
                )

            if not os.path.isfile(settings.MOTIF_CLUSTER_INFO):
                errors.append(
                    Warning('Motif cluster info not found.', id='querytgdb.W004')
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
