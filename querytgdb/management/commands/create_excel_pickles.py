#!/usr/bin/env python3
"""
This script reads a TargetDB output dataframe which is stored as a pickle object.
Also format the excel sheet, column sizes, colors etc.
"""

from django.core.management.base import BaseCommand

from ...utils.excel import create_excel_zip


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('-p', '--pickledir', help='Pickle Directory', required=False)

    def handle(self, *args, **options):
        create_excel_zip(options['pickledir'])
