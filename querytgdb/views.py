import gzip
import os
import shutil
import time
from functools import partial
from io import BytesIO
from threading import Lock

import matplotlib
from django.conf import settings
from django.core.files.storage import FileSystemStorage
from django.http import FileResponse, Http404, HttpResponseBadRequest, HttpResponseNotFound, JsonResponse
from django.utils.datastructures import MultiValueDictKeyError
from django.utils.decorators import method_decorator
from django.views.decorators.csrf import csrf_exempt
from django.views.generic import View
from pyparsing import ParseException

from querytgdb.utils.excel import create_export_zip
from .utils import GzipFileResponse, NetworkJSONEncoder, PandasJSONEncoder, cache_result, cache_view, convert_float, \
    metadata_to_dict, read_from_cache, svg_font_adder
from .utils.analysis_enrichment import AnalysisEnrichmentError, analysis_enrichment
from .utils.file import get_file, get_gene_lists, get_genes, get_network, merge_network_lists, network_to_filter_tfs, \
    network_to_lists
from .utils.formatter import format_data
from .utils.parser import QueryError, get_query_result
from .utils.summary import get_summary

# matplotlib import order issues
matplotlib.use('SVG')

import matplotlib.pyplot as plt

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ["DejaVu Sans"]

from querytgdb.utils.gene_list_enrichment import gene_list_enrichment, gene_list_enrichment_json
from .utils.motif_enrichment import NoEnrichedMotif, get_motif_enrichment_heatmap, get_motif_enrichment_json, \
    get_motif_enrichment_heatmap_table
from .utils.network import get_auc, get_network_json, get_network_stats

lock = Lock()

static_storage = FileSystemStorage(settings.QUERY_CACHE)
common_genes_storage = FileSystemStorage(settings.GENE_LISTS)


@method_decorator(csrf_exempt, name='dispatch')
class QueryView(View):
    def get(self, request, request_id):
        try:
            output = static_storage.path(f'{request_id}_pickle')

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

    def post(self, request, *args, **kwargs):
        try:
            request_id = request.POST['requestId']

            output = static_storage.path(f'{request_id}_pickle')
            os.makedirs(output, exist_ok=True)

            file_opts = {}

            targetgenes_file = get_file(request, "targetgenes", common_genes_storage)
            filter_tfs_file = get_file(request, "filtertfs")
            target_networks = get_file(request, "targetnetworks", common_genes_storage)

            if targetgenes_file:
                user_lists = get_gene_lists(targetgenes_file)

                if not target_networks:
                    cache_result(user_lists, f'{output}/target_genes.pickle.gz')  # cache the user list here

                file_opts["user_lists"] = user_lists

            if filter_tfs_file:
                file_opts["tf_filter_list"] = get_genes(filter_tfs_file)

            if target_networks:
                network = get_network(target_networks)

                try:
                    user_lists = file_opts["user_lists"]
                    # merge network with current user_lists
                    user_lists = merge_network_lists(user_lists, network)
                except KeyError:
                    user_lists = network_to_lists(network)

                file_opts["user_lists"] = user_lists

                graph_filter_list = network_to_filter_tfs(network)

                try:
                    tf_filter_list = file_opts["tf_filter_list"]
                    tf_filter_list = tf_filter_list.append(graph_filter_list)
                except KeyError:
                    tf_filter_list = graph_filter_list

                if not graph_filter_list.empty:
                    file_opts["tf_filter_list"] = tf_filter_list.drop_duplicates().sort_values()

                cache_result(network, f'{output}/target_network.pickle.gz')
                cache_result(user_lists, f'{output}/target_genes.pickle.gz')

            edges = request.POST.getlist('edges')

            # save the query
            with open(output + '/query.txt', 'w') as f:
                f.write(request.POST['query'].strip() + '\n')

            result, metadata, stats = get_query_result(query=request.POST['query'],
                                                       edges=edges,
                                                       cache_path=output,
                                                       size_limit=100000000,
                                                       **file_opts)

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
        except QueryError:
            return HttpResponseBadRequest("Result data size too large.")
        except ValueError as e:
            raise Http404('Query not available') from e
        except (MultiValueDictKeyError, ParseException):
            return HttpResponseBadRequest("Problem with query.")


class NetworkAuprView(View):
    def get(self, request, request_id):
        cache_dir = static_storage.path(f'{request_id}_pickle')

        precision = convert_float(request.GET.get('precision'))

        try:
            with lock:
                result = read_from_cache(get_auc)(
                    os.path.join(cache_dir, 'target_network.pickle.gz'),
                    os.path.join(cache_dir, 'tabular_output_unfiltered.pickle.gz'),
                    precision_cutoff=precision,
                    cache_path=cache_dir)

                return GzipFileResponse(result, content_type="image/svg+xml")
        except FileNotFoundError:
            raise Http404


class StatsView(View):
    def get(self, request, request_id):
        try:
            cache_dir = static_storage.path(f'{request_id}_pickle/tabular_output.pickle.gz')
            stats_cache = static_storage.path(f'{request_id}_pickle/stats.pickle.gz')

            info = cache_view(partial(read_from_cache(get_network_stats), cache_dir), stats_cache)

            return JsonResponse(info, encoder=PandasJSONEncoder)
        except FileNotFoundError:
            raise Http404


class NetworkJSONView(View):
    def get(self, request, request_id):
        try:
            edges = request.GET.getlist('edges')

            cache_dir = static_storage.path(f'{request_id}_pickle/tabular_output.pickle.gz')
            network_cache_dir = static_storage.path(f'{request_id}_pickle/network.pickle.gz')

            result = cache_view(
                partial(read_from_cache(get_network_json),
                        cache_dir),
                network_cache_dir, dummy_cache=True)

            return JsonResponse(result, safe=False, encoder=NetworkJSONEncoder)
        except ValueError:
            return HttpResponseBadRequest("Network too large", content_type="application/json")
        except FileNotFoundError:
            return HttpResponseNotFound(content_type="application/json")


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

            return FileResponse(open(out_file, 'rb'),
                                content_type='application/zip',
                                as_attachment=True,
                                filename='query.zip')

        except FileNotFoundError as e:
            return HttpResponseNotFound(content_type='application/zip')


class ListEnrichmentSVGView(View):
    def get(self, request, request_id):
        try:
            upper = convert_float(request.GET.get('upper'))
            lower = convert_float(request.GET.get('lower'))

            cache_path = static_storage.path("{}_pickle".format(request_id))

            buff = BytesIO()

            gene_list_enrichment(
                cache_path,
                draw=True,
                lower=lower,
                upper=upper
            ).savefig(buff)
            plt.close()

            buff.seek(0)
            svg_font_adder(buff)

            return FileResponse(buff, content_type='image/svg+xml')
        except (FileNotFoundError, ValueError) as e:
            return HttpResponseNotFound(content_type='image/svg+xml')


class ListEnrichmentLegendView(View):
    def get(self, request, request_id):
        try:
            cache_path = static_storage.path("{}_pickle".format(request_id))

            result = cache_view(
                partial(gene_list_enrichment, cache_path, legend=True),
                cache_path + '/list_enrichment_legend.pickle.gz'
            )
            return JsonResponse(result, safe=False, encoder=PandasJSONEncoder)
        except FileNotFoundError as e:
            raise Http404 from e


class ListEnrichmentTableView(View):
    def get(self, request, request_id):
        try:
            pickledir = static_storage.path("{}_pickle".format(request_id))

            result = cache_view(
                partial(gene_list_enrichment_json, pickledir),
                pickledir + '/list_enrichment.pickle.gz'
            )

            return JsonResponse(result, encoder=PandasJSONEncoder)
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

                cache_path = static_storage.path(f"{request_id}_pickle/tabular_output.pickle.gz")

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

                return FileResponse(buff, content_type='image/svg+xml')
            except (FileNotFoundError, NoEnrichedMotif):
                return HttpResponseNotFound(content_type='image/svg+xml')
            except (ValueError, TypeError, FloatingPointError):
                return HttpResponseBadRequest(content_type='image/svg+xml')


class MotifEnrichmentHeatmapTableView(View):
    def get(self, request, request_id):
        if not request_id:
            raise Http404

        cache_path = static_storage.path(f"{request_id}_pickle/tabular_output.pickle.gz")
        target_genes = static_storage.path(f"{request_id}_pickle/target_genes.pickle.gz")

        if not os.path.exists(cache_path):
            time.sleep(3)

        try:
            return JsonResponse(
                list(get_motif_enrichment_heatmap_table(
                    cache_path,
                    target_genes
                )),
                safe=False
            )
        except FileNotFoundError:
            raise Http404


class MotifEnrichmentInfo(View):
    def get(self, request):
        g = gzip.open(settings.MOTIF_CLUSTER)

        return GzipFileResponse(g,
                                content_type="text/csv",
                                filename="cluster_info.csv",
                                as_attachment=True)


class AnalysisEnrichmentView(View):
    def get(self, request, request_id):
        cache_path = static_storage.path(f"{request_id}_pickle/tabular_output.pickle.gz")
        analysis_cache = static_storage.path(f"{request_id}_pickle/analysis_enrichment.pickle.gz")

        try:
            result = cache_view(partial(analysis_enrichment, cache_path), analysis_cache)

            return JsonResponse(result, encoder=PandasJSONEncoder)
        except FileNotFoundError:
            return HttpResponseNotFound("Please make a new query")
        except AnalysisEnrichmentError as e:
            return HttpResponseBadRequest(e)


class SummaryView(View):
    def get(self, request, request_id):
        cache_path = static_storage.path(f"{request_id}_pickle/tabular_output.pickle.gz")
        summary_cache = static_storage.path(f"{request_id}_pickle/summary.pickle.gz")

        try:
            result = cache_view(partial(read_from_cache(get_summary), cache_path), summary_cache)

            return JsonResponse(result, encoder=PandasJSONEncoder)
        except FileNotFoundError:
            return HttpResponseNotFound("Please make a new query")
