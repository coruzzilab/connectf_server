import gzip
import mimetypes
import os
import os.path
import shutil

import yaml
from django.conf import settings
from django.core.management.base import BaseCommand, CommandError, CommandParser


def gz_copy(src, dst, force=False):
    mtype, enc = mimetypes.guess_type(src)

    if not force and os.path.isfile(dst):
        raise ValueError("Destination file exists, set 'force' to True to overwrite.")

    if mtype == 'text/csv':
        if enc is None:
            with open(src, 'rb') as f, gzip.open(dst, 'wb') as g:
                shutil.copyfileobj(f, g)

        elif enc == 'gzip':
            shutil.copyfile(src, dst)
    else:
        raise ValueError("Should be .csv file")


class Command(BaseCommand):
    def add_arguments(self, parser: CommandParser):
        parser.add_argument('-f', '--force', help="overwrite files even if they exist.", action='store_true')
        parser.add_argument("-d", "--description", help="motif description (.csv)", type=str)
        parser.add_argument("-m", "--motifs", help="gene to motif counts (.csv)", type=str)
        parser.add_argument("-i", "--individual-motifs", help="gene to individual motif counts (.csv)", type=str)

    def handle(self, *args, **options):
        if not (options['description'] or options['motifs'] or options['individual_motifs']):
            raise CommandError("Specify at least one of 'description', 'motifs', or 'individual-motifs'.")

        os.makedirs(os.path.join(settings.BASE_DIR, 'data'), exist_ok=True)  # make data directory if not exist

        opts = {}

        if options['description']:
            desc_path = os.path.join(settings.BASE_DIR, 'data/cluster_info.csv.gz')
            gz_copy(options['description'], desc_path, force=options['force'])
            opts['MOTIF_CLUSTER_INFO'] = desc_path

        if options['motifs']:
            motifs_path = os.path.join(settings.BASE_DIR, 'data/motifs.csv.gz')
            gz_copy(options['motifs'], motifs_path, force=options['force'])
            opts['MOTIF_ANNOTATION'] = motifs_path

        if options['individual_motifs']:
            motifs_indv_path = os.path.join(settings.BASE_DIR, 'data/motifs_indv.csv.gz')
            gz_copy(options['individual_motifs'], motifs_indv_path, force=options['force'])
            opts['MOTIF_TF_ANNOTATION'] = motifs_indv_path

        # Writes configs back into config.yaml with backup
        if opts:
            backup_path = settings.CONFIG_PATH + '.bak'
            shutil.copy2(settings.CONFIG_PATH, backup_path)
            self.stdout.write(f"Backup created at {backup_path}")
            new_configs = {**settings.CONFIG, **opts}
            with open(settings.CONFIG_PATH, 'w') as f:
                yaml.safe_dump(new_configs, f)

        self.stdout.write("Remember to restart the server.\n")
