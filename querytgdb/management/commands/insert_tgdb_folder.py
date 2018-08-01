import glob
import re

import os.path as path
import pandas as pd
from django.core.management.base import BaseCommand
from ...utils.insert_tgdb import insert_data

gene_regex = re.compile(r'.+(at[\dmc]g\d{5})(?:_MB030217_|_nlp7_|_tga1_)?(\d+N)?.+', re.I)


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('folder', help='path to tgdb data folder')

    def handle(self, *args, **options):
        data = pd.DataFrame(glob.glob(path.join(options['folder'], 'FDR_0.025_2/*.csv')))
        data[1] = data[0].str.replace(gene_regex, r'\1_\2')

        meta = pd.DataFrame(glob.glob(path.join(options['folder'], 'metadata/*.txt')))
        meta[1] = meta[0].str.replace(gene_regex, r'\1_\2')

        df = data.merge(meta, on=1, how='outer').sort_values(['0_y', '0_x'])

        for idx, row in df.iterrows():
            insert_data(row['0_x'], row['0_y'], None)
