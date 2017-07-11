#!/usr/bin/python
# -*- coding: utf-8 -*-

from ...models import TargetDBTF, Edges, Metadata, Analysis, Annotation, ReferenceId, \
    Interactions, Regulation, MetaIddata, AnalysisIddata, DAPdata
from django.core.management.base import BaseCommand
from django.db.models import Max
from decimal import Decimal
from collections import defaultdict
import sys

class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('-i', '--datafile', help= 'Input file for gene list')
        parser.add_argument('-m', '--metadata', help= 'Input file for metadata of an experiment')
        parser.add_argument('-d', '--dapdatafile', help= 'Input file for dap-seq data')


    def handle(self, *args, **options):
        print(options['datafile'])
        self.insertdata(options['datafile'], options['metadata'], options['dapdatafile'])

    # Parent function insertdata
    def insertdata(self, datafile, metafile, dapdatafile):

        # read target gene input file- store data in a list
        datalist = list()
        data = open(datafile, 'r')
        for vald in data:
            datalist.append(vald)

        # read meta data input file- store data in a dictionary
        metadict = dict()
        metadata = open(metafile, 'r')
        for valm in metadata:
            if not valm.split(':')[1].strip().upper() == 'NA':
                metadict[valm.split(':')[0].strip().upper()] = valm.split(':')[1].upper().strip()

        # Check if TF already exist in the database, if not then flag_dap=1;
        # flag_dap=1 means it is a new TF and DAP-seq data for this TF does not exist in db
        # call the function insert_dapseq to insert the dap-seq data for this TF
        flag_dap= 0
        alltfs = TargetDBTF.objects.values_list('db_tf_agi', flat=True)
        if not metadict['Transcription_Factor_ID'.upper()] in alltfs:
            flag_dap= 1

        #print('flag_dap= ',flag_dap)

        # Function to insert metadata
        curr_metaid, curr_analysisid, curr_refid= self.insert_metadata_analysis(metadict)

        # Functions to insert Nodes and Edges
        if metadict['Experiment_Type'.upper()] == 'Expression'.upper():
            if metadict['Experiment'.upper()] == 'Target'.upper():
                self.insert_nodes_target(datalist, metadict, curr_metaid, curr_analysisid, curr_refid)
        if metadict['Experiment_Type'.upper()] == 'Binding'.upper():
                self.insert_nodes_binding(datalist, metadict, curr_metaid, curr_analysisid, curr_refid)

        if flag_dap== 1:
            self.insert_dapseq(metadict['Transcription_Factor_ID'.upper()], dapdatafile)

    # FUNC1: insert metadata
    def insert_metadata_analysis(self, metadict):

        try:
            rs_meta= Metadata.objects.values_list('meta_fullid', flat=True) # get the existing meta_fullid from the database
            allmetaid= list(set(rs_meta))
        except Exception as e:
            allmetaid=[]

        meta_exist_flag = 0
        existed_exp_id = None

        # if experiment ID already exists: check if for the same exp, analysis_id provided in the current file also exists
        if metadict['Experiment_ID'.upper()].upper().strip() in allmetaid:
            print('Metaid already Exist!')
            meta_exist_flag = meta_exist_flag + 1
            ref_data= ReferenceId.objects.filter(meta_id__meta_fullid__exact=metadict['Experiment_ID'.upper()].upper(),
                                                analysis_id__analysis_fullid__exact=metadict['Analysis_ID'.upper()].upper())\
                                                .values('analysis_id')

            if not ref_data:
                exist_mdata_tmp = ReferenceId.objects.filter(meta_id__meta_fullid__exact=metadict['Experiment_ID'.upper()].
                                  upper()).distinct().values_list('meta_id', flat=True)
                existed_exp_id= list(exist_mdata_tmp)[0]
            else: # Cannot submit an experiment as both exp and analysis id already exist
                print('\nError: Experiment ID and analysis ID provided in the Metadata file ', metadict['Experiment_ID'.upper()].
                      upper(), ' and ', metadict['Analysis_ID'.upper()].upper(), ' already exists in the database')
                print('Provide a different Experiment ID or analysis ID to submit your data\nEXIT\n')
                sys.exit(1)

        # Getting the current metaid and analysisid
        meta_curr_val = metadict['Experiment_ID'.upper()]
        analysis_curr_val = metadict['Analysis_ID'.upper()]

        try:
            metaid_max = int(Metadata.objects.all().aggregate(Max('meta_id'))['meta_id__max'])
            # print('metaid_max= ',metaid_max)
            # max analysis id for a particular experiment
            analysisid_max = int(Analysis.objects.all().aggregate(Max('analysis_id'))['analysis_id__max'])
            refid_max = int(ReferenceId.objects.all().aggregate(Max('ref_id'))['ref_id__max'])
        except Exception as e:
            metaid_max = 0
            analysisid_max = 0
            refid_max = 0

        analysis_datatypes = ['ANALYSIS_ID', 'ANALYSIS_METHOD', 'ANALYSIS_CUTOFF', 'ANALYSIS_COMMAND', 'ANALYSIS NOTES']

        # inserts data to Metadata table:Experiment_ID as metaid
        if meta_exist_flag == 0:
            curr_metaid = metaid_max + 1
            curr_analysisid = analysisid_max + 1
            curr_refid = refid_max + 1
            mdata_obj = Metadata(meta_id=curr_metaid, meta_fullid=meta_curr_val)
            mdata_obj.save()
            adata_obj = Analysis(analysis_id=curr_analysisid, analysis_fullid=analysis_curr_val)
            adata_obj.save()
            for i in metadict:
                #print('i= ',i)
                if not i.upper() in analysis_datatypes:
                    #print('Not in analysis_datatypes')
                    middata_obj = MetaIddata(meta_id=mdata_obj, meta_value=metadict[i], meta_type=i)
                    middata_obj.save()
                else:
                    analysisdata_obj= AnalysisIddata(analysis_value=metadict[i], analysis_type=i, analysis_id= adata_obj)
                    analysisdata_obj.save()
            rdata_obj = ReferenceId(ref_id=curr_refid, meta_id=mdata_obj, analysis_id=adata_obj)
            rdata_obj.save()
        else:
            curr_metaid = existed_exp_id
            curr_analysisid = analysisid_max + 1
            curr_refid = refid_max + 1
            prev_mdata_obj= Metadata.objects.get(meta_id__exact=curr_metaid)
            adata_obj = Analysis(analysis_id=curr_analysisid, analysis_fullid=analysis_curr_val)
            adata_obj.save()
            for i in metadict:
                if i.upper() in analysis_datatypes:
                    analysisdata_obj = AnalysisIddata(analysis_value=metadict[i], analysis_type=i, analysis_id=adata_obj)
                    analysisdata_obj.save()
            rdata_obj = ReferenceId(ref_id=curr_refid, meta_id=prev_mdata_obj, analysis_id=adata_obj)
            rdata_obj.save()

        return curr_metaid, curr_analysisid, curr_refid


    # FUNC2: Target and inplanta- RNA-seq and microarray
    def insert_nodes_target(self, datalist, metadict, curr_metaid, curr_analysisid, curr_refid):

        # insert TF in TargetDBTf table, if TF already does not exist
        db_TF = TargetDBTF.objects.values_list('db_tf_agi', flat=True)

        #Insert TargetDB TF in TargetDBTF table
        exp_TF_name = metadict['Transcription_Factor_ID'.upper()]
        if not exp_TF_name in list(db_TF):
            tf_obj= TargetDBTF(db_tf_agi=exp_TF_name)
            tf_obj.save()

        edgelist = list()

        for i in datalist:
            edge_name = ':'.join([metadict['Experiment'.upper()], metadict['Expression_Type'.upper()],
                                  i.split('\t')[2].upper().strip(), i.split('\t')[3].upper().strip()])
            edgelist.append(edge_name)

        edgelist = list(set(edgelist))
        self.insert_edges(edgelist)

        TF_id = metadict['Transcription_Factor_ID'.upper()]
        db_tfid = list(TargetDBTF.objects.filter(db_tf_agi=TF_id).values_list('db_tf_id', flat= True))[0]
        #print('TF_id= ', TF_id, ', db_tfid= ', db_tfid)
        for k in datalist:
            tf_tg_edge = ':'.join([metadict['Experiment'.upper()], metadict['Expression_Type'.upper()],
                               k.split('\t')[2].upper().strip(), k.split('\t')[3].upper().strip()])
            inter_obj= Interactions(db_tf_id= TargetDBTF.objects.get(db_tf_agi__exact=TF_id),
                         target_id= Annotation.objects.get(agi_id= k.split('\t')[0]),
                         edge_id= Edges.objects.get(edge_name=tf_tg_edge),
                         ref_id= ReferenceId.objects.get(ref_id=curr_refid))
            inter_obj.save()

            pval = '%.1E' % Decimal(float(k.split('\t')[1].upper().strip()))
            fc = "%.1f" % float(k.split('\t')[4].upper().strip())

            reg_obj= Regulation(ref_id= ReferenceId.objects.get(ref_id=curr_refid),
                                ath_id= Annotation.objects.get(agi_id= k.split('\t')[0]),
                                foldchange=str(fc), pvalue=str(pval))
            reg_obj.save()


    # FUNC3: CHIP-SEQ and DAMID
    def insert_nodes_binding(self, datalist, metadict, curr_metaid, curr_analysisid, curr_refid):

        # insert TF in TargetDBTf table, if TF already does not exist
        db_TF = TargetDBTF.objects.values_list('db_tf_agi', flat=True)

        #Insert TargetDB TF in TargetDBTF table
        exp_TF_name = metadict['Transcription_Factor_ID'.upper()]
        if not exp_TF_name in list(db_TF):
            tf_obj= TargetDBTF(db_tf_agi=exp_TF_name)
            tf_obj.save()

        edgelist = list()

        for i in datalist:
            edge_name = ':'.join([metadict['Experiment'.upper()], metadict['Binding_Type'.upper()],
                                  i.split('\t')[2].upper().strip(), i.split('\t')[3].upper().strip()])
            edgelist.append(edge_name)

        edgelist = list(set(edgelist))
        self.insert_edges(edgelist)

        TF_id = metadict['Transcription_Factor_ID'.upper()]
        db_tfid = list(TargetDBTF.objects.filter(db_tf_agi=TF_id).values_list('db_tf_id', flat= True))[0]
        #print('TF_id= ', TF_id, ', db_tfid= ', db_tfid)
        for k in datalist:
            tf_tg_edge = ':'.join([metadict['Experiment'.upper()], metadict['Binding_Type'.upper()],
                               k.split('\t')[2].upper().strip(), k.split('\t')[3].upper().strip()])
            inter_obj= Interactions(db_tf_id= TargetDBTF.objects.get(db_tf_agi__exact=TF_id),
                         target_id= Annotation.objects.get(agi_id= k.split('\t')[0]),
                         edge_id= Edges.objects.get(edge_name=tf_tg_edge),
                         ref_id= ReferenceId.objects.get(ref_id=curr_refid))
            inter_obj.save()

            #pval = '%.1E' % Decimal(float(k.split('\t')[1].upper().strip()))
            #fc = "%.1f" % float(k.split('\t')[4].upper().strip())

            #reg_obj= Regulation(ref_id= ReferenceId.objects.get(ref_id=curr_refid),
            #                    ath_id= Annotation.objects.get(agi_id= k.split('\t')[0]),
            #                    foldchange=str(fc), pvalue=str(pval))
            #reg_obj.save()


    # FUNC4: Insert edges in Edges table
    def insert_edges(self, edgelist):
        all_edges = Edges.objects.values_list('edge_name', flat= True)
        for j in edgelist:
            if not j in list(all_edges):
                edge_obj= Edges(edge_name=j)
                edge_obj.save()


    # FUNC5: Insert DAP-seq data if TF is already not in the database
    def insert_dapseq(self, tfid, dapdatafile):

        dapdata = open(dapdatafile, 'r')

        # get the gene number from the Annotation table
        daptf = list()
        daptargets= list()
        dapinteract = defaultdict(list)
        for valt in dapdata:
            if valt.split('\t')[0].strip() == tfid:
                daptf.append(valt.split('\t')[0].strip())
                daptargets.append(valt.split('\t')[1].strip())
                dapinteract[valt.split('\t')[0].strip()].append(valt.split('\t')[1].strip())

        #print('length= ', dapinteract)

        daptf_number= {z[0]: z[1] for z in
                           list(TargetDBTF.objects.filter(db_tf_agi__in=daptf).values_list('db_tf_agi', 'db_tf_id'))}

        daptargets_number = {x[0]: x[1] for x in
                           list(Annotation.objects.filter(agi_id__in=daptargets).values_list('agi_id', 'ath_id'))}
        #print('daptargets_number= ', daptargets_number)

        for val_d in dapinteract:
            #print('val_d= ',val_d)
            #print('daptf_number[val_d]= ',daptf_number[val_d])
            for val_d_tg in dapinteract[val_d]:

                dapdata_obj = DAPdata(db_tfid=TargetDBTF.objects.get(db_tf_id__exact=
                                                                  daptf_number[val_d]),
                                      ath_id=Annotation.objects.get(ath_id__exact=
                                                                  daptargets_number[val_d_tg]))

                dapdata_obj.save()


