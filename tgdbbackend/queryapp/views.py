from django.shortcuts import render
import os;
from django.views.generic import View;
from django.views.decorators.csrf import csrf_exempt
from django.utils.decorators import method_decorator
from django.http import JsonResponse
from .app import query_TargetDB;
import tempfile;
import shutil;
import pandas as pd;
# Create your views here.
def save_file(dest_path, f):
    #original_name, file_extension = os.path.splitext(f.name)
    #filename = filename + '-' + datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S') + file_extension
    path = os.path.join( dest_path , f.name);
    destination = open(path, 'wb+')
    for chunk in f.chunks():
        destination.write(chunk)
    destination.close()
    return path

@method_decorator(csrf_exempt, name='dispatch')
class HandleQueryView(View):
    
    def get(self,request,*args,**kwargs):
        pass
    

    def post(self,request,*args,**kwargs):
        request_id = request.POST['requestId'];
        dbname = 'targetdb';
        TFquery = request.POST['tfs'].split(" ") if request.POST['tfs']!='' else None;
        edges = request.POST['edges'].split(" ") if request.POST['edges']!='' else None;
        metadata = request.POST['metas'].split(" ") if request.POST['metas']!='' else None;
        
        targetgenesFilePath = None;
        dirpath = tempfile.mkdtemp();
        
        tfFilePaths = [];
        if len(request.FILES)!=0:
            if request.FILES.has_key("targetgenes"):
                targetgenesFilePath = save_file(dirpath, request.FILES["targetgenes"]) 
            if request.FILES.has_key("file-0"):
                i=0;
                while request.FILES.has_key("file-"+str(i)):
                    tfFilePaths.append(save_file(dirpath, request.FILES["file-"+str(i)]))
                    i+=1
        
        self.setTFquery(TFquery,tfFilePaths);
        
        outfmt = 'tabular';       
        output = None;
        #df,out_metadata_df = query_TargetDB.main(dbname, TFquery, edges, metadata, outfmt, output, 
                           #targetgenesFilePath);
        df = pd.read_pickle("testdf.pickle")
        out_metadata_df = pd.read_pickle("testmeta.pickle")     
        
        df_columns = [];
        for column in df.columns:
            df_columns.append({"id": column, "name": column, "field": column});
        
        res = [{'columns':df_columns}];
        res[0]['data'] = df.to_json(orient='index');
        meta_columns=[];
        out_metadata_df.reset_index(inplace= True);
        for column in out_metadata_df.columns:
            meta_columns.append({"id": column, "name": column, "field": column});
        
        res.append({'columns':meta_columns, 'data':out_metadata_df.to_json(orient='index')})
        shutil.rmtree(dirpath);
        return JsonResponse(res,safe=False);
    
    #----------------------------------------------------------------------
    def setTFquery(self,TFquery,tfFilePaths):
        """"""
        i = 0;
        j = 0;
        if (TFquery==None): return;
        for i in range(len(TFquery)):
            if TFquery[i].find("{")!=-1:
                begin  = TFquery[i].find("{")
                end = TFquery[i].find("}")
                TFquery[i] = TFquery[i][0:begin]+ tfFilePaths[j]+TFquery[i][end+1:]
                j+=1