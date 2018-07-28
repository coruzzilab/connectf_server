import os
import pathlib
import re
import shutil
import time
from io import BytesIO, TextIOWrapper
from threading import Lock
from typing import List

import pandas as pd
from django.conf import settings
from django.core.files.storage import FileSystemStorage, SuspiciousFileOperation
from django.http import Http404, HttpResponse, HttpResponseBadRequest, HttpResponseNotFound, JsonResponse
from django.utils.datastructures import MultiValueDictKeyError
from django.utils.decorators import method_decorator
from django.views.decorators.csrf import csrf_exempt
from django.views.generic import View
from pyparsing import ParseException

from querytgdb.models import ReferenceId
from querytgdb.utils.clustering import heatmap
from querytgdb.utils.excel import create_export_zip
from .utils import CytoscapeJSONEncoder, PandasJSONEncoder, cache_result, convert_float, metadata_to_dict
from .utils.cytoscape import get_cytoscape_json
from .utils.file import get_gene_lists
from .utils.formatter import format_data
from .utils.motif_enrichment import NoEnrichedMotif, get_motif_enrichment_heatmap, get_motif_enrichment_json
from .utils.parser import get_query_result

lock = Lock()

APPS_DIR = pathlib.Path(settings.BASE_DIR) / 'tgdbbackend'
STATIC_DIR = APPS_DIR / 'static' / 'queryBuilder'

static_storage = FileSystemStorage(pathlib.Path(settings.MEDIA_ROOT) / 'queryBuilder')
common_genes_storage = FileSystemStorage('commongenelists')


def parse_heatmap_idx(idx: str) -> List[str]:
    return re.split(r'[\s|]+', idx)[1:-1]


@method_decorator(csrf_exempt, name='dispatch')
class QueryView(View):
    def post(self, request, *args, **kwargs):
        try:
            request_id = request.POST['requestId']

            output = static_storage.path(request_id + '_pickle')
            os.makedirs(output, exist_ok=True)

            targetgenes_file = None

            if 'targetgenes' in request.POST:
                try:
                    targetgenes_file = common_genes_storage.open("{}.txt".format(request.POST['targetgenes']), 'r')
                except (FileNotFoundError, SuspiciousFileOperation):
                    pass

            if request.FILES:
                if "targetgenes" in request.FILES:
                    targetgenes_file = TextIOWrapper(request.FILES["targetgenes"])

            if targetgenes_file:
                user_lists = get_gene_lists(targetgenes_file)
                result, metadata, stats = get_query_result(request.POST['query'],
                                                           user_lists=user_lists,
                                                           cache_path=output)
            else:
                result, metadata, stats = get_query_result(request.POST['query'], cache_path=output)

            columns, merged_cells, result_list = format_data(result, stats)

            cache_result(result_list, output + '/formatted_tabular_output.pickle.gz')

            res = [
                {
                    'data': result_list,
                    'mergeCells': merged_cells,
                    'columns': columns
                },
                metadata_to_dict(metadata)
            ]

            return JsonResponse(res, safe=False, encoder=PandasJSONEncoder)
        except ValueError as e:
            raise Http404('Query not available') from e
        except (MultiValueDictKeyError, ParseException):
            return HttpResponseBadRequest("Propblem with query.")


class QueryIdView(View):
    def get(self, request, request_id):
        try:
            output = static_storage.path(request_id + '_pickle')

            result, metadata, stats = get_query_result(cache_path=output)

            columns, merged_cells, result_list = format_data(result, stats)

            res = [
                {
                    'data': result_list,
                    'mergeCells': merged_cells,
                    'columns': columns
                },
                metadata_to_dict(metadata)
            ]

            return JsonResponse(res, safe=False, encoder=PandasJSONEncoder)
        except FileNotFoundError as e:
            raise Http404('Query not available') from e


class StatsView(View):
    def get(self, request, request_id):
        try:
            cache_dir = static_storage.path(request_id + '_pickle')
            df = pd.read_pickle(cache_dir + '/tabular_output.pickle.gz')

            return JsonResponse({
                'num_edges': df.loc[:, (slice(None), slice(None), slice(None), 'EDGE')].count().sum(),
                'num_targets': df.shape[0]
            }, encoder=PandasJSONEncoder)
        except FileNotFoundError:
            raise Http404


class CytoscapeJSONView(View):
    def get(self, request, request_id):
        try:
            cache_dir = static_storage.path(request_id + '_pickle')
            df = pd.read_pickle(cache_dir + '/tabular_output.pickle.gz')
            return JsonResponse(get_cytoscape_json(df), safe=False, encoder=CytoscapeJSONEncoder)
        except ValueError as e:
            return HttpResponseBadRequest("Network too large")
        except FileNotFoundError as e:
            raise Http404 from e


class FileExportView(View):
    def get(self, request, request_id):
        try:
            if not request_id:
                raise FileNotFoundError

            out_file = static_storage.path("{}.zip".format(request_id))
            if not os.path.exists(out_file):
                cache_folder = static_storage.path("{}_pickle".format(request_id))
                out_folder = static_storage.path(request_id)

                shutil.rmtree(out_folder, ignore_errors=True)
                os.makedirs(out_folder)

                create_export_zip(cache_folder, out_folder)

                shutil.make_archive(out_folder, 'zip', out_folder)  # create a zip file for output directory
                shutil.rmtree(out_folder, ignore_errors=True)  # delete the output directory after creating zip file

            with open(out_file, 'rb') as f:
                response = HttpResponse(f, content_type='application/zip')
                response['Content-Disposition'] = 'attachment; filename="query.zip"'

                return response

        except FileNotFoundError as e:
            return HttpResponseNotFound(content_type='application/zip')


class HeatMapPNGView(View):
    def get(self, request, request_id):
        try:
            upper = convert_float(request.GET.get('upper'))
            lower = convert_float(request.GET.get('lower'))

            cache_path = static_storage.path("{}_pickle".format(request_id))

            buff = BytesIO()
            response = HttpResponse(content_type='image/svg+xml')

            heatmap(
                cache_path,
                draw=True,
                lower=lower,
                upper=upper
            ).savefig(buff)

            buff.seek(0)
            shutil.copyfileobj(buff, response)

            return response
        except (FileNotFoundError, ValueError) as e:
            return HttpResponseNotFound(content_type='image/svg+xml')


class HeatMapTableView(View):
    def get(self, request, request_id):
        try:
            df = heatmap(
                static_storage.path("{}_pickle".format(request_id)),
                draw=False
            )
            meta_ids, analysis_ids = zip(*map(parse_heatmap_idx, df.index))
            references = ReferenceId.objects.filter(analysis_id__analysis_fullid__in=analysis_ids,
                                                    meta_id__meta_fullid__in=meta_ids)

            def get_rows():
                for (idx, *row), analysis_id, meta_id in zip(df.itertuples(name=None), analysis_ids, meta_ids):
                    info = {'name': idx}
                    try:
                        ref = references.get(
                            analysis_id__analysis_fullid=analysis_id,
                            meta_id__meta_fullid=meta_id)
                        info.update(
                            ref.analysis_id.analysisiddata_set.values_list(
                                'analysis_type',
                                'analysis_value'))
                        info.update(ref.meta_id.metaiddata_set.values_list('meta_type', 'meta_value'))
                    except ReferenceId.DoesNotExist:
                        pass

                    yield [info] + row

            return JsonResponse({
                'columns': df.columns,
                'result': list(get_rows())
            }, encoder=PandasJSONEncoder)
        except (FileNotFoundError, ValueError) as e:
            raise Http404 from e


class MotifEnrichmentView(View):
    def get(self, request, request_id):
        if not request_id:
            return JsonResponse({}, status=404)
        with lock:
            try:
                alpha = float(request.GET.get('alpha', 0.05))
                body = request.GET.get('body', '0')

                cache_path = static_storage.path("{}_pickle/tabular_output.pickle.gz".format(request_id))

                if not os.path.exists(cache_path):
                    time.sleep(3)

                return JsonResponse(
                    get_motif_enrichment_json(
                        cache_path,
                        # static_storage.path("{}_pickle/target_genes.pickle.gz".format(request_id)),
                        alpha=alpha,
                        body=body == '1'),
                    encoder=PandasJSONEncoder)
            except (FileNotFoundError, NoEnrichedMotif) as e:
                return JsonResponse({}, status=404)
            except (ValueError, TypeError) as e:
                return JsonResponse({'error': str(e)}, status=400)


class MotifEnrichmentHeatmapView(View):
    def get(self, request, request_id):
        if not request_id:
            return HttpResponseNotFound(content_type='image/svg+xml')
        with lock:
            try:
                alpha = float(request.GET.get('alpha', 0.05))
                body = request.GET.get('body', '0')

                cache_path = static_storage.path("{}_pickle/tabular_output.pickle.gz".format(request_id))

                if not os.path.exists(cache_path):
                    time.sleep(3)

                buff = get_motif_enrichment_heatmap(
                    cache_path,
                    # static_storage.path("{}_pickle/target_genes.pickle.gz".format(request_id)),
                    alpha=alpha,
                    body=body == '1'
                )

                response = HttpResponse(content_type='image/svg+xml')
                shutil.copyfileobj(buff, response)

                return response
            except (FileNotFoundError, NoEnrichedMotif):
                return HttpResponseNotFound(content_type='image/svg+xml')
            except (ValueError, TypeError, FloatingPointError):
                return HttpResponseBadRequest(content_type='image/svg+xml')
