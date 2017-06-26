import os
import shutil
import tempfile
from functools import partial
from itertools import chain, groupby

import environ
import pandas as pd
from django.http import Http404, HttpResponse, JsonResponse
from django.utils.decorators import method_decorator
from django.views.decorators.csrf import csrf_exempt
from django.views.generic import View

from querytgdb.utils import query_tgdb
from querytgdb.utils.cytoscape import create_cytoscape_data
from querytgdb.utils.excel import create_excel_zip
from querytgdb.utils.clustering import read_pickled_targetdbout
from .utils import PandasJSONEncoder

ROOT_DIR = environ.Path(
    __file__) - 3  # (tgdbbackend/config/settings/common.py - 3 = tgdbbackend/)
APPS_DIR = ROOT_DIR.path('tgdbbackend')
STATIC_DIR = APPS_DIR.path('static').path('queryBuilder')


# Create your views here.
def save_file(dest_path, f):
    # original_name, file_extension = os.path.splitext(f.name)
    # filename = filename + '-' + datetime.datetime.now().strftime(
    # '%Y-%m-%d-%H-%M-%S') + file_extension
    path = os.path.join(dest_path, f.name)
    destination = open(path, 'wb+')
    for chunk in f.chunks():
        destination.write(chunk)
    destination.close()
    return path


@method_decorator(csrf_exempt, name='dispatch')
class HandleQueryView(View):
    def get(self, request, *args, **kwargs):
        pass

    def post(self, request, *args, **kwargs):
        request_id = request.POST['requestId']
        # dbname = 'targetdb'
        tf_query = request.POST['tfs'].split(" ") if request.POST[
                                                         'tfs'] != '' else None
        edges = request.POST['edges'].split(" ") if request.POST[
                                                        'edges'] != '' else \
            None
        metadata = request.POST['metas'].split(" ") if request.POST[
                                                           'metas'] != '' \
            else None

        targetgenes_file_path = None
        dirpath = tempfile.mkdtemp()

        tf_file_paths = []
        if len(request.FILES) != 0:
            if "targetgenes" in request.FILES:
                targetgenes_file_path = save_file(dirpath,
                                                  request.FILES["targetgenes"])
            if "file-0" in request.FILES:
                i = 0
                while "file-" + str(i) in request.FILES:
                    tf_file_paths.append(
                        save_file(dirpath, request.FILES["file-" + str(i)]))
                    i += 1

        self.setTFquery(tf_query, tf_file_paths)

        output = STATIC_DIR.path(request_id)
        df, out_metadata_df = query_tgdb(tf_query, edges, metadata, targetgenes_file_path, str(output))
        # df = pd.read_pickle("testdf.pickle")
        # out_metadata_df = pd.read_pickle("testmeta.pickle")

        # df_columns = [{"id": column, "name": column, "field": column} for column in df.columns]

        # res = [{'columns': df_columns, 'data': df.to_json(orient='index')}]

        # print(df.columns)
        int_cols = df.columns.get_level_values(2).isin(['Pvalue', 'Foldchange']) | df.columns.get_level_values(0).isin(
            ['UserList', 'Target Count'])
        df.iloc[4:, int_cols] = df.iloc[4:, int_cols].apply(partial(pd.to_numeric, errors='coerce'))
        df['Ind'] = df['Ind'].astype(int)
        df = df.where(pd.notnull(df), None)

        merged_cells = []

        for i, level in enumerate(df.columns.labels):
            index = 0
            for label, group in groupby(level):
                size = sum(1 for _ in group)
                merged_cells.append({'row': i, 'col': index, 'colspan': size, 'rowspan': 1})
                if i == 0:
                    merged_cells.extend({'row': a, 'col': index, 'colspan': size, 'rowspan': 1} for a in range(3, 6))
                index += size

        res = [{
            'data': list(chain(zip(*df.columns), df.itertuples(index=False, name=None))),
            'mergeCells': merged_cells
        }]

        out_metadata_df.reset_index(inplace=True)
        meta_columns = [{"id": column, "name": column, "field": column} for column in out_metadata_df.columns]

        res.append({'columns': meta_columns,
                    'data': out_metadata_df.to_json(orient='index')})
        shutil.rmtree(dirpath)

        return JsonResponse(res, safe=False, encoder=PandasJSONEncoder)

    # ----------------------------------------------------------------------
    def setTFquery(self, TFquery, tfFilePaths):
        i = 0
        j = 0
        if TFquery is None:
            return
        for i in range(len(TFquery)):
            if TFquery[i].find("{") != -1:
                begin = TFquery[i].find("{")
                end = TFquery[i].find("}")
                TFquery[i] = TFquery[i][0:begin] + tfFilePaths[j] + TFquery[i][
                                                                    end + 1:]
                j += 1


class CytoscapeJSONView(View):
    def get(self, request, request_id, name):
        try:
            outdir = str(STATIC_DIR.path("{}_json".format(request_id)))
            if not os.path.isdir(outdir):
                outdir = create_cytoscape_data(str(STATIC_DIR.path("{}_pickle".format(request_id))))
            with open("{}/{}.json".format(outdir, name)) as f:
                return HttpResponse(f, content_type="application/json; charset=utf-8")
        except FileNotFoundError as e:
            raise Http404 from e


class ExcelDownloadView(View):
    def get(self, request, request_id):
        try:
            out_file = str(STATIC_DIR.path("{}.zip".format(request_id)))
            if not os.path.exists(out_file):
                out_file = create_excel_zip(str(STATIC_DIR.path("{}_pickle".format(request_id))))
            with open(out_file, 'rb') as f:
                response = HttpResponse(f, content_type='application/zip')
                response['Content-Disposition'] = 'attachment; filename="query.zip"'

                return response
        except FileNotFoundError as e:
            raise Http404 from e


class HeatMapPNGView(View):
    def get(self, request, request_id):
        try:
            response = HttpResponse(content_type='image/svg+xml')
            read_pickled_targetdbout(
                str(STATIC_DIR.path("{}_pickle".format(request_id))),
                save_file=False
            ).savefig(response)

            return response
        except FileNotFoundError as e:
            raise Http404 from e
