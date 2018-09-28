import os
import pathlib
import shutil
import time
from io import BytesIO, TextIOWrapper
from threading import Lock

import matplotlib
import pandas as pd
from django.conf import settings
from django.core.files.storage import FileSystemStorage, SuspiciousFileOperation
from django.http import Http404, HttpResponse, HttpResponseBadRequest, HttpResponseNotFound, JsonResponse
from django.utils.datastructures import MultiValueDictKeyError
from django.utils.decorators import method_decorator
from django.views.decorators.csrf import csrf_exempt
from django.views.generic import View
from pyparsing import ParseException

from querytgdb.models import Analysis
from querytgdb.utils.excel import create_export_zip
from .utils import CytoscapeJSONEncoder, PandasJSONEncoder, cache_result, convert_float, metadata_to_dict, \
    svg_font_adder
from .utils.analysis_enrichment import analysis_enrichment
from .utils.cytoscape import get_cytoscape_json
from .utils.file import get_gene_lists
from .utils.formatter import format_data
from .utils.parser import get_query_result
from .utils.summary import get_summary

# matplotlib import order issues
matplotlib.use('SVG')

import matplotlib.pyplot as plt

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ["DejaVu Sans"]

from querytgdb.utils.gene_list_enrichment import gene_list_enrichment
from .utils.motif_enrichment import NoEnrichedMotif, get_motif_enrichment_heatmap, get_motif_enrichment_json, \
    get_motif_enrichment_heatmap_table

lock = Lock()

APPS_DIR = pathlib.Path(settings.BASE_DIR) / 'tgdbbackend'
STATIC_DIR = APPS_DIR / 'static' / 'queryBuilder'

static_storage = FileSystemStorage(pathlib.Path(settings.MEDIA_ROOT) / 'queryBuilder')
common_genes_storage = FileSystemStorage('commongenelists')


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

            edges = request.POST.getlist('edges')

            if targetgenes_file:
                user_lists = get_gene_lists(targetgenes_file)
                result, metadata, stats = get_query_result(request.POST['query'],
                                                           user_lists=user_lists,
                                                           edges=edges,
                                                           cache_path=output)
            else:
                result, metadata, stats = get_query_result(request.POST['query'],
                                                           edges=edges,
                                                           cache_path=output)

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

            info = {
                'num_edges': df.loc[:, (slice(None), slice(None), slice(None), 'EDGE')].count().sum(),
                'num_targets': df.shape[0]
            }

            return JsonResponse(info, encoder=PandasJSONEncoder)
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


class ListEnrichmentSVGView(View):
    def get(self, request, request_id):
        try:
            upper = convert_float(request.GET.get('upper'))
            lower = convert_float(request.GET.get('lower'))

            cache_path = static_storage.path("{}_pickle".format(request_id))

            buff = BytesIO()
            response = HttpResponse(content_type='image/svg+xml')

            gene_list_enrichment(
                cache_path,
                draw=True,
                lower=lower,
                upper=upper
            ).savefig(buff)

            buff.seek(0)
            svg_font_adder(buff)
            buff.seek(0)
            shutil.copyfileobj(buff, response)

            return response
        except (FileNotFoundError, ValueError) as e:
            return HttpResponseNotFound(content_type='image/svg+xml')


class ListEnrichmentLegendView(View):
    def get(self, request, request_id):
        try:
            cache_path = static_storage.path("{}_pickle".format(request_id))
            return JsonResponse(gene_list_enrichment(
                cache_path,
                legend=True
            ), safe=False, encoder=PandasJSONEncoder)
        except FileNotFoundError as e:
            raise Http404 from e


class ListEnrichmentTableView(View):
    def get(self, request, request_id):
        try:
            df = gene_list_enrichment(
                static_storage.path("{}_pickle".format(request_id)),
                draw=False
            )
            names, criteria, exp_ids, analysis_ids, ls, uids = zip(*df.index)
            analyses = Analysis.objects.filter(name__in=analysis_ids, experiment__name__in=exp_ids)

            def get_rows():
                for (name, criterion, exp_id, analysis_id, l, uid), *row in df.itertuples(name=None):
                    info = {'name': name, 'filter': criterion, 'targets': l}
                    try:
                        analysis = analyses.get(
                            name=analysis_id,
                            experiment__name=exp_id)
                        info.update(
                            analysis.analysisdata_set.values_list('key', 'value'))
                        info.update(analysis.experiment.experimentdata_set.values_list('key', 'value'))
                    except Analysis.DoesNotExist:
                        pass

                    yield [info] + row

            return JsonResponse({
                'columns': df.columns,
                'result': list(get_rows())
            }, encoder=PandasJSONEncoder)
        except (FileNotFoundError, ValueError) as e:
            raise Http404 from e


class MotifEnrichmentJSONView(View):
    def get(self, request, request_id):
        if not request_id:
            raise Http404
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
                        static_storage.path("{}_pickle/target_genes.pickle.gz".format(request_id)),
                        alpha=alpha,
                        body=body == '1'),
                    encoder=PandasJSONEncoder)
            except (FileNotFoundError, NoEnrichedMotif) as e:
                raise Http404 from e
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
                upper = convert_float(request.GET.get('upper'))
                lower = convert_float(request.GET.get('lower'))

                cache_path = static_storage.path("{}_pickle/tabular_output.pickle.gz".format(request_id))

                if not os.path.exists(cache_path):
                    time.sleep(3)

                buff = get_motif_enrichment_heatmap(
                    cache_path,
                    static_storage.path("{}_pickle/target_genes.pickle.gz".format(request_id)),
                    upper_bound=upper,
                    lower_bound=lower,
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


class MotifEnrichmentHeatmapTableView(View):
    def get(self, request, request_id):
        if not request_id:
            raise Http404

        cache_path = static_storage.path("{}_pickle/tabular_output.pickle.gz".format(request_id))

        if not os.path.exists(cache_path):
            time.sleep(3)

        try:
            return JsonResponse(
                list(get_motif_enrichment_heatmap_table(
                    cache_path,
                    static_storage.path("{}_pickle/target_genes.pickle.gz".format(request_id))
                )),
                safe=False
            )
        except FileNotFoundError:
            raise Http404


class AnalysisEnrichmentView(View):
    def get(self, request, request_id):
        cache_path = static_storage.path("{}_pickle/tabular_output.pickle.gz".format(request_id))

        try:
            return JsonResponse(analysis_enrichment(cache_path), encoder=PandasJSONEncoder)

        except FileNotFoundError:
            return HttpResponseNotFound("Please make a new query")
        except ValueError as e:
            return HttpResponseBadRequest(e)


class SummaryView(View):
    def get(self, request, request_id):
        cache_path = static_storage.path("{}_pickle/tabular_output.pickle.gz".format(request_id))

        try:
            result = get_summary(cache_path)

            return JsonResponse(result, encoder=PandasJSONEncoder)
        except FileNotFoundError:
            return HttpResponseNotFound("Please make a new query")
