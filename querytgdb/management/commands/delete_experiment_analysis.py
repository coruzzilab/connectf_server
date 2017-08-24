#!/usr/bin/env python

'''
This script is to discard records (experiments or analysis) from the targetdb database.
-e option: If only experiment id is provided, it will delete the whole experiment (means all the analysis in an exp).
-a option: If analysis id is provided along with the experiment id, it will delete the individual analysis.


'''


import os,sys
from ...models import Metadata, Analysis, ReferenceId, Interactions, Regulation
from django.core.management.base import BaseCommand


class Command(BaseCommand):

    def add_arguments(self, parser):
        # program only deletes one experiment/analysis at a time, so the argument provided are not lists
        parser.add_argument('-e', '--experimentid', help= 'delete an experiment', required= True)
        parser.add_argument('-a', '--analysisid', help= 'delete an analysis within an experiment')

    def handle(self, *args, **options):
        self.delete_data(options['experimentid'], options['analysisid'])

    def delete_data(self, userexperimentid, useranalysisid):
        if userexperimentid and useranalysisid:  # delete an analysis
            self.del_analysis(userexperimentid, useranalysisid)

        if not useranalysisid:  # delete the whole experiment
            self.del_exp(userexperimentid)


    # delete an analysis
    def del_analysis(self, userexperimentid, useranalysisid):
        print('userexperimentid= ',userexperimentid)
        print('useranalysisid= ',useranalysisid)

        # Get the reference id for the experimentid and analysisid
        rs_refid = Metadata.objects.select_related('referenceid').filter(meta_fullid__exact = userexperimentid,
                                    referenceid__analysis_id__analysis_fullid__exact = useranalysisid).\
                                    values_list('referenceid__ref_id','referenceid__meta_id','referenceid__analysis_id')


        print('rs_refid= ',rs_refid)

        # In case of deleting an analysis, only a particular reference id will be deleted

        # Data to be deleted from:
        # referenceid table, metadata, metaiddata, analysis, analysisiddata, regulation, interactions
        #ReferenceId.objects.filter(ref_id__in = rs_refid).delete()
        #Interactions.objects.filter(ref_id__in = rs_refid).delete()
        #Regulation.objects.filter(ref_id__in = rs_refid).delete()



    # delete an experiment
    def del_exp(self, userexperimentid):
        print('userexperimentid= ',userexperimentid)

        # Get the referenceid for user provided experimentid
        rs_refid = Metadata.objects.select_related('referenceid').filter(meta_fullid__exact = userexperimentid).\
                                     values_list('referenceid__ref_id','referenceid__meta_id','referenceid__analysis_id')


        print('rs_refid= ',rs_refid)

        # In case of deleting an experiment, all the reference ids belonging to an experiment will be deleted
        #ReferenceId.objects.filter(ref_id__in = rs_refid).delete()
        #Interactions.objects.filter(ref_id__in = rs_refid).delete()
        #Regulation.objects.filter(ref_id__in = rs_refid).delete()






