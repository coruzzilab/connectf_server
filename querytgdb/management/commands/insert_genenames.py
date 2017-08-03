#!/usr/bin/python
# -*- coding: utf-8 -*-

''' This script inserts the target DB data into the database

'''


##############
# Modules
##############

import sys
from django.core.management.base import BaseCommand
from django.db.models import Max
from ...models import Annotation

class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('datafile')

    def handle(self, *args, **options):
        print(options['datafile'])
        self.insert_data(options['datafile'])

    def insert_data(self, datafile):

        tairdata= open(datafile,'r', encoding='utf-8')

        for valt in tairdata:
            anno_obj= Annotation(agi_id=valt.split('\t')[0].upper().strip(),
                                 ath_gene_type=valt.split('\t')[1].upper().strip(),
                                 ath_name=valt.split('\t')[2].upper().strip(),
                                 ath_fullname=valt.split('\t')[3].upper().strip(),
                                 ath_gene_fam=valt.split('\t')[4].upper().strip())

            anno_obj.save()
