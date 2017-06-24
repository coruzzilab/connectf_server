'''
This module returns four json objects:

    Query TF- Query TF
	Query TF- database TF
	database TF- whole genome TF
	database TF- whole genome (TF+nonTF)

Notes: create_json is the master function. Returns three objects.

Last updated: April 12, 2017
'''

##############
# Modules

from __future__ import absolute_import

from django.core.management.base import BaseCommand

from ...utils.cytoscape import create_cytoscape_data


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('-p', '--pickledir', help='Pickle Directory', required=False)

    def handle(self, *args, **options):
        create_cytoscape_data(options['pickledir'])
