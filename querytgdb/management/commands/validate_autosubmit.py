"""
This script validates the data based on the crietria below:
    1) If metadata file entires are correct
    2) All the uploaded files are in tab-delimited format
    3) Have correct number of columns and data entered in files is in specified format
    4) No duplicated genes in read count and gene list file
    5) All the samples provided in read count should be included in experimental design file
If the data pass all the above validation steps, run the insert script.
"""
from django.core.management.base import BaseCommand

from ...utils.validate_autosubmit import validate_genelist, validate_metadata, validate_readcount_expdesign


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('-i', '--genelist', help='User submitted genelist file', required=True)
        parser.add_argument('-m', '--metadatafile', help='User submitted metadata file', required=True)
        parser.add_argument('-r', '--readcount', help='User submitted readcount file', required=True)
        parser.add_argument('-e', '--expdesign', help='User submitted experimental design file', required=True)

    def handle(self, *args, **options):
        metadict = validate_metadata(options['metadatafile'])
        validate_genelist(options['genelist'], metadict)
        validate_readcount_expdesign(options['readcount'], options['expdesign'])
