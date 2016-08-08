from django.shortcuts import render

from django.views.generic import View;

# Create your views here.
class HandleQueryView(View):
    
    def get(self,request,*args,**kwargs):
        print (1);
        pass
    
    def post(self,request,*args,**kwargs):
        pass;