from django.shortcuts import render

from django.views.generic import View;
from django.views.decorators.csrf import csrf_exempt
from django.utils.decorators import method_decorator
from django.http import JsonResponse
from .app import query_TargetDB;

import pandas as pd;
# Create your views here.
@method_decorator(csrf_exempt, name='dispatch')
class HandleQueryView(View):
    
    def get(self,request,*args,**kwargs):
        pass
    

    def post(self,request,*args,**kwargs):
        dbname = 'targetdb';
        TFquery = ['AT4G24020'];
        TFquery = request.POST['tfs'].split(" ") if request.POST['tfs']!='' else None;
        edges = request.POST['edges'].split(" ") if request.POST['edges']!='' else None;
        metadata = request.POST['metas'].split(" ") if request.POST['metas']!='' else None;
        outfmt = 'tabular';
        targetgenes = request.POST['targetgenes'] if request.POST['targetgenes']!='' else None;
        output = None;
        #df,out_metadata_df = query_TargetDB.main(dbname, TFquery, edges, metadata, outfmt, output, 
                           #targetgenes);
        df = pd.read_pickle("testdf.pickle")
        out_metadata_df = pd.read_pickle("testmeta.pickle")        
        # data.columns is the list of column
        # Index([u'Gene Full Name', u'Gene Name', u'Gene ID', u'Target Count', u'p-value', u'AT2G22200_MB041416_RNASEQ', u'AT4G24020_AS050515_CHIPSEQ', u'AT4G24020_AS081114_RNASEQ', u'AT4G24020_AS090116_RNASEQ'], dtype='object')
        df_columns = [];
        for column in df.columns:
            df_columns.append({"id": column, "name": column, "field": column});
        
        # data.index is the rows
        # new_res_df[0:5].to_json(orient='index')
        # '{"0":{"Gene Full Name":null,"Gene Name":null,"Gene ID":"Targets","Target Count":null,"p-value":null,"AT2G22200_MB041416_RNASEQ":"DREB (2209)","AT4G24020_AS050515_CHIPSEQ":"NLP7 (5693)","AT4G24020_AS081114_RNASEQ":"NLP7 (3692)","AT4G24020_AS090116_RNASEQ":"NLP7 (751)"},"1":{"Gene Full Name":null,"Gene Name":null,"Gene ID":"-","Target Count":null,"p-value":null,"AT2G22200_MB041416_RNASEQ":null,"AT4G24020_AS050515_CHIPSEQ":"0:5:10:30:180","AT4G24020_AS081114_RNASEQ":null,"AT4G24020_AS090116_RNASEQ":null},"2":{"Gene Full Name":"MYB DOMAIN PROTEIN 74","Gene Name":"MYB74,ATMYB74","Gene ID":"AT4G05100","Target Count":4.0,"p-value":null,"AT2G22200_MB041416_RNASEQ":"INDUCED:PCHX","AT4G24020_AS050515_CHIPSEQ":"11111","AT4G24020_AS081114_RNASEQ":"REPRESSED:MCHX","AT4G24020_AS090116_RNASEQ":"REPRESSED:PCHX"}}'
        res = [{'columns':df_columns}];
        res[0]['data'] = df.to_json(orient='index');
        meta_columns=[];
        for column in out_metadata_df.columns:
            meta_columns.append({"id": column, "name": column, "field": column});
        res.append({'columns':meta_columns, 'data':out_metadata_df.to_json(orient='index')})
        return JsonResponse(res,safe=False);