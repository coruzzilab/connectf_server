#!/usr/bin/python
# -*- coding: utf-8 -*-

from django.core.management.base import BaseCommand, CommandError

from ...utils import query_tgdb


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('-t', '--TFname', nargs='+', help='Search by TF name or get all the data '
                                                              'from the database (-t= OR [ALLTF]', required=True)
        parser.add_argument('-e', '--edges', nargs='+', help='search by edges')
        parser.add_argument('-m', '--metadata', nargs='+', help='search by metadata')
        parser.add_argument('-o', '--output', help='output file name', required=False)
        parser.add_argument('-r', '--targetgenes',
                            help='list of genes provided by the user to refine the database output')


    def handle(self, *args, **options):
        try:
            query_tgdb(options['TFname'], options['edges'], options['metadata'],
                       options['targetgenes'], options['output'])
        except TypeError as e:
            raise CommandError(str(e)) from e
