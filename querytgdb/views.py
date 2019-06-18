import logging
import os
import shutil
import tempfile
import warnings
from threading import Lock
from uuid import uuid4

import matplotlib
from django.conf import settings
from django.core.cache import cache
from django.core.files.storage import FileSystemStorage
from django.http import FileResponse, Http404, HttpResponse, HttpResponseBadRequest, HttpResponseNotFound, JsonResponse
from django.utils.datastructures import MultiValueDictKeyError
from django.views.generic import View

from querytgdb.utils.export import create_export_zip, export_csv, write_excel
from querytgdb.utils.gene_list_enrichment import gene_list_enrichment, gene_list_enrichment_json
from .utils import GzipFileResponse, NetworkJSONEncoder, PandasJSONEncoder, check_annotations, convert_float, \
    metadata_to_dict, svg_font_adder
from .utils.analysis_enrichment import AnalysisEnrichmentError, analysis_enrichment
from .utils.file import BadFile, get_file, get_gene_lists, get_genes, get_network, merge_network_lists, \
    network_to_filter_tfs, network_to_lists
from .utils.formatter import format_data
from .utils.motif_enrichment import MOTIF, NoEnrichedMotif, get_motif_enrichment_heatmap, \
    get_motif_enrichment_heatmap_table, get_motif_enrichment_json
from .utils.network import get_auc_figure, get_network_json, get_network_sif, get_network_stats, get_pruned_network
from .utils.parser import QueryError, get_query_result
from .utils.summary import get_summary

logger = logging.getLogger(__name__)

lock = Lock()

gene_lists_storage = FileSystemStorage(settings.GENE_LISTS)
networks_storage = FileSystemStorage(settings.TARGET_NETWORKS)


class QueryView(View):
    def get(self, request, request_id):
        cached_data = cache.get(f'{request_id}/formatted_tabular_output')
        if cached_data is None:
            raise Http404('Query not available')

        columns, merged_cells, result_list, metadata = cached_data

        res = {
            'result': {
                'data': result_list,
                'mergeCells': merged_cells,
                'columns': columns,
            },
            'metadata': metadata,
            'request_id': request_id
        }

        return JsonResponse(res, encoder=PandasJSONEncoder)

    def post(self, request, *args, **kwargs):
        errors = []

        try:
            request_id = str(uuid4())

            file_opts = {}

            targetgenes_file = get_file(request, "targetgenes", gene_lists_storage)
            filter_tfs_file = get_file(request, "filtertfs")
            target_networks = get_file(request, "targetnetworks", networks_storage)

            if targetgenes_file:
                user_lists = get_gene_lists(targetgenes_file)

                bad_genes = check_annotations(user_lists[0].index)
                if bad_genes:
                    errors.append(f'Genes in Target Genes File not in database: {", ".join(bad_genes)}')

                if not target_networks:
                    cache.set(f'{request_id}/target_genes', user_lists)  # cache the user list here

                file_opts["user_lists"] = user_lists

            if filter_tfs_file:
                filter_tfs = get_genes(filter_tfs_file)
                file_opts["tf_filter_list"] = filter_tfs

                bad_genes = check_annotations(filter_tfs)
                if bad_genes:
                    errors.append(f'Genes in Filter TFs File not in database: {", ".join(bad_genes)}')

            if target_networks:
                network = get_network(target_networks)

                bad_genes = check_annotations(network[1]['source'].append(network[1]['target']))
                if bad_genes:
                    errors.append(f'Genes in Network File not in database: {", ".join(bad_genes)}')

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

                cache.set_many({f'{request_id}/target_network': network,
                                f'{request_id}/target_genes': user_lists})

            edges = request.POST.getlist('edges')
            query = request.POST['query']

            # save the query
            cache.set(f'{request_id}/query', query.strip() + '\n')

            result, metadata, stats, _uid = get_query_result(query=query,
                                                             edges=edges,
                                                             size_limit=50_000_000,
                                                             uid=request_id,
                                                             **file_opts)

            metadata = metadata_to_dict(metadata)
            columns, merged_cells, result_list = format_data(result, stats)

            cache.set(f'{request_id}/formatted_tabular_output', (columns, merged_cells, result_list, metadata))

            res = {
                'result': {
                    'data': result_list,
                    'mergeCells': merged_cells,
                    'columns': columns
                },
                'metadata': metadata,
                'request_id': request_id
            }

            if errors:
                res['errors'] = errors

            return JsonResponse(res, encoder=PandasJSONEncoder)
        except (QueryError, BadFile) as e:
            return HttpResponseBadRequest(e)
        except ValueError as e:
            raise Http404(f'Query not available: {e}') from e
        except MultiValueDictKeyError:
            return HttpResponseBadRequest("Problem with query.")


class NetworkAuprView(View):
    def get(self, request, request_id):
        precision = convert_float(request.GET.get('precision'))

        try:
            with lock, warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=matplotlib.MatplotlibDeprecationWarning,
                                        module="matplotlib.figure")

                cached_data = cache.get_many([f'{request_id}/target_network',
                                              f'{request_id}/tabular_output_unfiltered'])

                result = get_auc_figure(cached_data[f'{request_id}/target_network'],
                                        cached_data[f'{request_id}/tabular_output_unfiltered'],
                                        request_id,
                                        precision_cutoff=precision)

                return GzipFileResponse(result, content_type="image/svg+xml")
        except KeyError:
            raise Http404

    def head(self, request, request_id):
        if cache.get(f'{request_id}/target_network') is not None:
            return HttpResponse()
        raise Http404


class NetworkPrunedView(View):
    def get(self, request, request_id, cutoff):
        try:
            try:
                precision_cutoff = float(cutoff)
            except (ValueError, TypeError):
                return HttpResponseBadRequest("Invalid precision cutoff")

            pruned = get_pruned_network(request_id, precision_cutoff)

            response = HttpResponse(content_type='text/csv')
            response['Content-Disposition'] = 'attachment; filename="pruned_network.csv"'

            pruned.to_csv(response, index=False)

            return response
        except ValueError:
            return HttpResponseNotFound(content_type='text/csv')


class StatsView(View):
    def get(self, request, request_id):
        try:
            data_cache = f'{request_id}/tabular_output'
            stats_cache = f'{request_id}/stats'

            stats = cache.get(stats_cache)

            if stats is None:
                try:
                    stats = get_network_stats(cache.get(data_cache))
                except AttributeError:
                    raise Http404

                cache.set(stats_cache, stats)

            return JsonResponse(stats, encoder=PandasJSONEncoder)
        except FileNotFoundError:
            raise Http404


class NetworkJSONView(View):
    def get(self, request, request_id):
        try:
            edges = request.GET.getlist('edges')
            precision = convert_float(request.GET.get('precision'))

            result = get_network_json(request_id, edges=edges, precision_cutoff=precision)

            return JsonResponse(result, safe=False, encoder=NetworkJSONEncoder)
        except ValueError:
            return HttpResponseBadRequest("Network too large", content_type="application/json")
        except KeyError:
            return HttpResponseNotFound(content_type="application/json")


class NetworkSifView(View):
    def get(self, request, request_id):
        try:
            edges = request.GET.getlist('edges')
            precision = convert_float(request.GET.get('precision'))

            resp = HttpResponse(content_type='text/plain')
            resp['Content-Disposition'] = f'attachment; filename="{request_id}.sif"'

            get_network_sif(request_id, edges=edges, precision_cutoff=precision, buffer=resp)

            return resp
        except ValueError:
            return HttpResponseBadRequest("Network too large", content_type="text/plain")
        except KeyError:
            return HttpResponseNotFound(content_type="text/plain")


class FileExportView(View):
    def get(self, request, request_id):
        try:
            out_file = os.path.join(tempfile.gettempdir(), f'{request_id}.zip')
            if not os.path.exists(out_file):

                with tempfile.TemporaryDirectory() as temp_folder:
                    create_export_zip(request_id, temp_folder)
                    shutil.make_archive(os.path.splitext(out_file)[0], 'zip', temp_folder,
                                        logger=logger)  # create a zip file for output directory

            return FileResponse(open(out_file, 'rb'),
                                content_type='application/zip',
                                as_attachment=True,
                                filename='query.zip')

        except KeyError:
            return HttpResponseNotFound(content_type='application/zip')


class ExcelExportView(View):
    def get(self, request, request_id):
        try:
            with tempfile.NamedTemporaryFile(suffix='.xlsx') as f:
                write_excel(request_id, f.name)
                f.seek(0)

                response = HttpResponse(
                    content_type='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet')
                response['Content-Disposition'] = 'attachment; filename="output.xlsx"'
                shutil.copyfileobj(f, response)

                return response
        except KeyError:
            return HttpResponseNotFound(
                content_type='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet')


class CsvExportView(View):
    def get(self, request, request_id):
        response = HttpResponse(content_type='text/csv')
        response['Content-Disposition'] = 'attachment; filename="output.csv"'

        try:
            df = export_csv(request_id)
            df.to_csv(response)

            return response
        except KeyError:
            return HttpResponseNotFound(content_type='text/csv')


class ListEnrichmentSVGView(View):
    def get(self, request, request_id):
        try:
            upper = convert_float(request.GET.get('upper'))
            lower = convert_float(request.GET.get('lower'))

            buff = gene_list_enrichment(
                request_id,
                draw=True,
                lower=lower,
                upper=upper
            )
            svg_font_adder(buff)

            return FileResponse(buff, content_type='image/svg+xml')
        except ValueError:
            return HttpResponseNotFound(content_type='image/svg+xml')


class ListEnrichmentLegendView(View):
    def get(self, request, request_id):
        try:
            result = cache.get(f'{request_id}/list_enrichment_legend')

            if result is None:
                result = gene_list_enrichment(request_id, legend=True)
                cache.set(f'{request_id}/list_enrichment_legend', result)

            return JsonResponse(result, safe=False, encoder=PandasJSONEncoder)
        except ValueError as e:
            raise Http404 from e


class ListEnrichmentTableView(View):
    def get(self, request, request_id):
        try:
            result = cache.get(f'{request_id}/list_enrichment')

            if result is None:
                result = gene_list_enrichment_json(request_id)
                cache.set(f'{request_id}/list_enrichment', result)

            return JsonResponse(result, encoder=PandasJSONEncoder)
        except ValueError as e:
            raise Http404(e) from e


class MotifEnrichmentJSONView(View):
    def get(self, request, request_id):
        if not request_id:
            raise Http404
        with lock:
            try:
                try:
                    alpha = float(request.GET.get('alpha', 0.05))
                except ValueError:
                    alpha = 0.05
                regions = request.GET.getlist('regions')

                return JsonResponse(
                    get_motif_enrichment_json(
                        request_id,
                        regions,
                        alpha=alpha),
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
                try:
                    alpha = float(request.GET.get('alpha', 0.05))
                except ValueError:
                    alpha = 0.05
                regions = request.GET.getlist('regions')
                upper = convert_float(request.GET.get('upper'))
                lower = convert_float(request.GET.get('lower'))

                buff = get_motif_enrichment_heatmap(
                    request_id,
                    regions,
                    upper_bound=upper,
                    lower_bound=lower,
                    alpha=alpha
                )

                return FileResponse(buff, content_type='image/svg+xml')
            except (FileNotFoundError, NoEnrichedMotif):
                return HttpResponseNotFound(content_type='image/svg+xml')
            except (ValueError, TypeError, FloatingPointError):
                return HttpResponseBadRequest(content_type='image/svg+xml')


class MotifEnrichmentHeatmapTableView(View):
    def get(self, request, request_id):
        try:
            return JsonResponse(
                list(get_motif_enrichment_heatmap_table(
                    request_id
                )),
                safe=False
            )
        except KeyError:
            raise Http404


class MotifEnrichmentInfo(View):
    def get(self, request):
        return GzipFileResponse(open(settings.MOTIF_CLUSTER, 'rb'),
                                content_type="text/csv",
                                filename="cluster_info.csv",
                                as_attachment=True)


class MotifEnrichmentRegions(View):
    def get(self, request):
        return JsonResponse({
            'regions': MOTIF.region_desc,
            'default_regions': MOTIF.default_regions
        })


class AnalysisEnrichmentView(View):
    def get(self, request, request_id):
        try:
            result = cache.get(f'{request_id}/analysis_enrichment')

            if result is None:
                df = cache.get(f'{request_id}/tabular_output')

                if df is None:
                    return HttpResponseNotFound("Please make a new query")

                result = analysis_enrichment(df)
                cache.set(f'{request_id}/analysis_enrichment', result)

            return JsonResponse(result, encoder=PandasJSONEncoder)
        except AnalysisEnrichmentError as e:
            return HttpResponseBadRequest(e)


class SummaryView(View):
    def get(self, request, request_id):
        result = cache.get(f'{request_id}/summary')

        if result is None:
            df = cache.get(f"{request_id}/tabular_output")
            if df is None:
                return HttpResponseNotFound("Please make a new query")

            result = get_summary(df)
            cache.set(f'{request_id}/summary', result)

        return JsonResponse(result, encoder=PandasJSONEncoder)
