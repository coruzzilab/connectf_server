#!/usr/bin/python
# -*- coding: utf-8 -*-

from django.core.management.base import BaseCommand

from ...utils.insert_tgdb import insertdata


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('-i', '--datafile', help='Input file for gene list')
        parser.add_argument('-m', '--metadata', help='Input file for metadata of an experiment')
        parser.add_argument('-d', '--dapdatafile', help='Input file for dap-seq data')

    def handle(self, *args, **options):
        insertdata(options['datafile'], options['metadata'], options['dapdatafile'])
