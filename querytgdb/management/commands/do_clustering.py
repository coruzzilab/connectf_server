#!/usr/bin/env python3

"""
This script 1) reads TargetDB output dataframe and targetgene list that are stored as a pickle object
            2) Intersects the dataframe list with the targetgenelist
            3) creates a matrix of number of overlaps between different lists
            4) perform enrichment-test on the matrix and creates a new matrix with p-values
            5) converts matrix p-values to heatmaps
"""

from django.core.management.base import BaseCommand
from ...utils.clustering import heatmap


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('-p', '--pickledir', help='Pickle Directory', required=False)
        parser.add_argument('-c', '--cutoff', help='cutoff number of genes for enrichment test, default= 10', required=False)
        parser.add_argument('-b', '--background', help='background Araport11/TAIR10, default= 28775', required=False)

    def handle(self, *args, **options):
        heatmap(options['pickledir'], options['cutoff'], options['background'])
