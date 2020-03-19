import json
import logging
import mimetypes
import os
import re
import shutil
import tempfile
import warnings
from itertools import chain
from operator import itemgetter
from threading import Lock
from typing import Optional
from uuid import uuid4

import matplotlib
from django.conf import settings
from django.core.cache import cache
from django.core.files.storage import FileSystemStorage
from django.http import FileResponse, Http404, HttpResponse, HttpResponseBadRequest, HttpResponseNotFound, JsonResponse
from django.utils.datastructures import MultiValueDictKeyError
from django.views.generic import View
from jsonschema import ValidationError, validate

from querytgdb.utils.export import create_export_zip, export_csv, write_excel
from querytgdb.utils.gene_list_enrichment import gene_list_enrichment, gene_list_enrichment_json
from .utils import GzipFileResponse, NetworkJSONEncoder, PandasJSONEncoder, check_annotations, \
    convert_float, metadata_to_dict, svg_font_adder
from .utils.analysis_enrichment import AnalysisEnrichmentError, analysis_enrichment, analysis_enrichment_csv
from .utils.file import BadFile, filter_gene_lists_by_background, get_background_genes, get_file, get_gene_lists, \
    get_genes, get_network, merge_network_filter_tfs, merge_network_lists, network_to_filter_tfs, network_to_lists
from .utils.formatter import format_data
from .utils.motif_enrichment import ADD_MOTIFS, MOTIFS, MotifEnrichmentError, NoEnrichedMotif, \
    get_additional_motif_enrichment_json, get_motif_enrichment_heatmap, get_motif_enrichment_heatmap_table, \
    get_motif_enrichment_json
from .utils.motif_enrichment.motif import AdditionalMotifData, MotifData
from .utils.network import get_auc_figure, get_network_json, get_network_sif, get_network_stats, get_pruned_network
from .utils.parser import Ids, QueryError, filter_df_by_ids, get_query_result, reorder_data
from .utils.summary import get_summary

logger = logging.getLogger(__name__)

lock = Lock()  # Ensure all matplotlib functions use this lock

gene_lists_storage = FileSystemStorage(settings.GENE_LISTS)
networks_storage = FileSystemStorage(settings.TARGET_NETWORKS)


class QueryView(View):
    """
    Enpoint for new query or get cached queries
    """

    def get(self, request, request_id):
        try:
            cached_data = cache.get_many([f'{request_id}/formatted_tabular_output', f'{request_id}/analysis_ids'])

            columns, merged_cells, result_list, metadata = cached_data[f'{request_id}/formatted_tabular_output']
            ids = cached_data[f'{request_id}/analysis_ids']

            res = {
                'result': {
                    'data': result_list,
                    'mergeCells': merged_cells,
                    'columns': columns,
                },
                'metadata': metadata,
                'request_id': request_id,
                'analysis_ids': list(ids.items())
            }

            return JsonResponse(res, encoder=PandasJSONEncoder)
        except KeyError:
            try:
                result, metadata, stats, _uid, ids = get_query_result(size_limit=50_000_000,
                                                                      uid=request_id)

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
                    'request_id': request_id,
                    'analysis_ids': list(ids.items())
                }

                return JsonResponse(res, encoder=PandasJSONEncoder)
            except KeyError:
                raise Http404('Query not available')

    def post(self, request, *args, **kwargs):
        errors = []

        try:
            request_id = str(uuid4())

            file_opts = {}

            targetgenes_file, targetgenes_source = get_file(request, "targetgenes", gene_lists_storage)
            filter_tfs_file, filter_tfs_source = get_file(request, "filtertfs")
            target_networks, networks_source = get_file(request, "targetnetworks", networks_storage)
            background_genes_file, background_genes_source = get_file(request, "backgroundgenes")

            if background_genes_file:
                background_genes = get_background_genes(background_genes_file)
                file_opts['target_filter_list'] = background_genes
                cache.set(f'{request_id}/background_genes', background_genes)

            if targetgenes_file:
                user_lists = get_gene_lists(targetgenes_file)

                if 'target_filter_list' in file_opts:
                    user_lists = filter_gene_lists_by_background(user_lists, file_opts['target_filter_list'])

                bad_genes = check_annotations(user_lists[0].index)
                if bad_genes and targetgenes_source != 'storage':
                    errors.append(f'Genes in Target Genes File not in database: {", ".join(bad_genes)}')

                if not target_networks:
                    cache.set(f'{request_id}/target_genes', user_lists)  # cache the user list here

                file_opts["user_lists"] = user_lists

            if filter_tfs_file:
                filter_tfs = get_genes(filter_tfs_file)
                file_opts["tf_filter_list"] = filter_tfs

                bad_genes = check_annotations(filter_tfs)
                if bad_genes and filter_tfs_source != 'storage':
                    errors.append(f'Genes in Filter TFs File not in database: {", ".join(bad_genes)}')

            if target_networks:
                network = get_network(target_networks)

                bad_genes = check_annotations(network[1]['source'].append(network[1]['target']))
                if bad_genes and networks_source != 'storage':
                    errors.append(f'Genes in Network File not in database: {", ".join(bad_genes)}')

                try:
                    user_lists = file_opts["user_lists"]
                    # merge network with current user_lists
                    network, user_lists = merge_network_lists(network, user_lists)
                except KeyError:
                    user_lists = network_to_lists(network)

                file_opts["user_lists"] = user_lists

                try:
                    tf_filter_list = file_opts["tf_filter_list"]
                    network, tf_filter_list = merge_network_filter_tfs(network, tf_filter_list)
                except KeyError:
                    tf_filter_list = network_to_filter_tfs(network)

                file_opts["tf_filter_list"] = tf_filter_list

                cache.set_many({f'{request_id}/target_network': network,
                                f'{request_id}/target_genes': user_lists})

            edges = request.POST.getlist('edges')
            query = request.POST['query']

            result, metadata, stats, _uid, ids = get_query_result(query=query,
                                                                  edges=edges,
                                                                  size_limit=50_000_000,
                                                                  uid=request_id,
                                                                  **file_opts)

            metadata = metadata_to_dict(metadata)
            columns, merged_cells, result_list = format_data(result, stats)

            cache.set_many({
                f'{request_id}/query': query.strip() + '\n',  # save queries
                f'{request_id}/formatted_tabular_output': (columns, merged_cells, result_list, metadata),
            })

            res = {
                'result': {
                    'data': result_list,
                    'mergeCells': merged_cells,
                    'columns': columns
                },
                'metadata': metadata,
                'request_id': request_id,
                'analysis_ids': list(ids.items())
            }

            if errors:
                res['errors'] = errors

            return JsonResponse(res, encoder=PandasJSONEncoder)
        except (QueryError, BadFile) as e:
            return HttpResponseBadRequest(e)
        except ValueError as e:
            return HttpResponseNotFound(f'Query not available: {e}', content_type='text/plain')
        except MultiValueDictKeyError as e:
            return HttpResponseBadRequest(f"Problem with query: {e}")


ANALYSIS_ID_SCHEMA = {
    "type": "array",
    "items": {
        "type": "array",
        "items": [
            {
                "type": "array",
                "items": [
                    {"type": "string"},
                    {"type": "number"}
                ]
            },
            {
                "type": "object",
                "properties": {
                    "show": {"type": "boolean"},
                    "name": {"type": "string"}
                },
                "additionalProperties": False,
                "required": ["show", "name"]
            }
        ]
    }
}


class EditQueryView(View):
    """
    Endpoints to hide or rename queried analyses
    """

    def get(self, request, request_id):
        ids = cache.get(f'{request_id}/analysis_ids')

        if ids is None:
            return HttpResponseNotFound('Query Analysis Ids not available')

        return JsonResponse(list(ids.items()), safe=False, encoder=PandasJSONEncoder)

    def post(self, request, request_id):
        try:
            data = json.loads(request.body)
            validate(data, ANALYSIS_ID_SCHEMA)

            ids: Optional[Ids] = cache.get(f'{request_id}/analysis_ids')

            if ids is None:
                return HttpResponseNotFound('Query Analysis Ids not available')

            data = {tuple(idx): opt for idx, opt in data}

            if not any(v['show'] for v in data.values()):
                return HttpResponseBadRequest("cannot hide all analyses")

            for key in ids:
                try:
                    ids[key] = data[key]
                except KeyError:
                    pass

            cache.set(f'{request_id}/analysis_ids', ids)

            cached_result = cache.get_many([f'{request_id}/target_genes', f'{request_id}/tabular_output_unfiltered'])

            result = cached_result[f'{request_id}/tabular_output_unfiltered']
            result = filter_df_by_ids(result, ids)

            try:
                user_lists = cached_result[f'{request_id}/target_genes']
                result = result[result.index.str.upper().isin(user_lists[0].index.str.upper())].dropna(axis=1,
                                                                                                       how='all')

                if result.empty:
                    raise QueryError("Empty result (user list too restrictive).")

                result = reorder_data(result)
            except KeyError:
                pass

            cache.set(f'{request_id}/tabular_output', result)  # refresh filtered tabular output

            # delete cache keys and refresh cache here.
            cache.delete_many([
                # formatted output
                f'{request_id}/formatted_tabular_output',
                # network
                f'{request_id}/network',
                f'{request_id}/network_data',
                # AUPR curve
                f'{request_id}/figure',
                f'{request_id}/figure_data',
                # network stats
                f'{request_id}/stats',
                # Gene list enrichment
                f'{request_id}/list_enrichment',
                f'{request_id}/list_enrichment_legend',
                # motif enrichment
                *(f'{request_id}/{r}_enrich' for r in MOTIFS.regions),
                # analysis enrichment
                f'{request_id}/analysis_enrichment',
                # summary
                f'{request_id}/summary'
            ])

            return JsonResponse(list(ids.items()), status=201, safe=False, encoder=PandasJSONEncoder)
        except (json.JSONDecodeError, ValidationError, QueryError) as e:
            return HttpResponseBadRequest(e)
        except KeyError:
            return HttpResponseNotFound('Query does not exist. Please start a new query.')


class NetworkAuprView(View):
    def get(self, request, request_id):
        precision = convert_float(request.GET.get('precision'))

        try:
            with lock, warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=matplotlib.MatplotlibDeprecationWarning,
                                        module="matplotlib.figure")

                cached_data = cache.get_many([
                    f'{request_id}/target_network',
                    f'{request_id}/tabular_output_unfiltered',
                    f'{request_id}/analysis_ids'
                ])

                df, ids = itemgetter(
                    f'{request_id}/tabular_output_unfiltered',
                    f'{request_id}/analysis_ids'
                )(cached_data)
                df = filter_df_by_ids(df, ids)

                result = get_auc_figure(cached_data[f'{request_id}/target_network'],
                                        df,
                                        request_id,
                                        precision_cutoff=precision)

                return GzipFileResponse(result, content_type="image/svg+xml")
        except KeyError:
            raise Http404

    def head(self, request, request_id):
        if cache.get(f'{request_id}/target_network') is not None:
            return HttpResponse()
        return HttpResponseNotFound()


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
                    raise Http404('Stats not available')

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
            expand = request.GET.get('expand')

            resp = HttpResponse(content_type='text/plain')
            resp['Content-Disposition'] = f'attachment; filename="{request_id}.sif"'

            get_network_sif(request_id, edges=edges, precision_cutoff=precision, expand=expand, buffer=resp)

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
            export_csv(request_id, response)

            return response
        except KeyError:
            return HttpResponseNotFound(content_type='text/csv')


class ListEnrichmentHeatmapView(View):
    def get(self, request, request_id):
        try:
            with lock:
                upper = convert_float(request.GET.get('upper'))
                lower = convert_float(request.GET.get('lower'))
                fields = request.GET.getlist('fields')

                buff = gene_list_enrichment(
                    request_id,
                    draw=True,
                    lower=lower,
                    upper=upper,
                    fields=fields
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

                background_genes = cache.get(f'{request_id}/background_genes')

                if background_genes is not None:
                    motif_data = MotifData(background=background_genes)
                else:
                    motif_data = MOTIFS

                return JsonResponse(
                    get_motif_enrichment_json(
                        request_id,
                        regions,
                        alpha=alpha,
                        motif_data=motif_data),
                    encoder=PandasJSONEncoder)
            except (FileNotFoundError, NoEnrichedMotif, KeyError) as e:
                raise Http404 from e
            except (MotifEnrichmentError, ValueError, TypeError) as e:
                return JsonResponse({'error': str(e)}, status=400)


class AdditionalMotifEnrichmentJSONView(View):
    def post(self, request, request_id):
        try:
            opts = {}
            data = json.loads(request.body)
            regions = data.get('regions', [])
            if 'motifs' in data:
                opts['motifs'] = data['motifs']

            background_genes = cache.get(f'{request_id}/background_genes')

            if background_genes is not None:
                motif_data = AdditionalMotifData(background=background_genes)
            else:
                motif_data = ADD_MOTIFS

            return JsonResponse(
                get_additional_motif_enrichment_json(
                    request_id,
                    regions,
                    use_default_motifs=True,
                    motif_data=motif_data,
                    **opts),
                encoder=PandasJSONEncoder)
        except (FileNotFoundError, NoEnrichedMotif, KeyError) as e:
            return JsonResponse({'error': str(e)}, status=404)
        except (MotifEnrichmentError, ValueError, TypeError, json.JSONDecodeError) as e:
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
                fields = request.GET.getlist('fields')

                background_genes = cache.get(f'{request_id}/background_genes')

                if background_genes is not None:
                    motif_data = MotifData(background=background_genes)
                else:
                    motif_data = MOTIFS

                buff = get_motif_enrichment_heatmap(
                    request_id,
                    regions,
                    upper_bound=upper,
                    lower_bound=lower,
                    alpha=alpha,
                    fields=fields,
                    motif_data=motif_data
                )

                return FileResponse(buff, content_type='image/svg+xml')
            except (FileNotFoundError, NoEnrichedMotif, KeyError):
                return HttpResponseNotFound(content_type='image/svg+xml')
            except (MotifEnrichmentError, ValueError, TypeError, FloatingPointError):
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
        with open(settings.MOTIF_CLUSTER, 'rb') as f:
            response = HttpResponse(content_type="text/csv")
            response["Content-Disposition"] = 'attachment; filename="cluster_info.csv"'
            response["Content-Encoding"] = 'gzip'

            shutil.copyfileobj(f, response)

            return response


class MotifEnrichmentRegions(View):
    def get(self, request):
        return JsonResponse({
            'regions': MOTIFS.region_desc,
            'default_regions': MOTIFS.default_regions
        })


class MotifEnrichmentMotifs(View):
    def get(self, request):
        return JsonResponse({
            'motifs': MOTIFS.motifs
        }, encoder=PandasJSONEncoder)


class MotifEnrichmentAdditionalMotifs(View):
    def get(self, request):
        return JsonResponse({
            'motifs': ADD_MOTIFS.motifs
        }, encoder=PandasJSONEncoder)


class AnalysisEnrichmentView(View):
    def get(self, request, request_id):
        try:
            result = cache.get(f'{request_id}/analysis_enrichment')

            if result is None:
                result = analysis_enrichment(request_id)
                cache.set(f'{request_id}/analysis_enrichment', result)

            return JsonResponse(result, encoder=PandasJSONEncoder)
        except AnalysisEnrichmentError as e:
            return HttpResponseBadRequest(e)


class AnalysisEnrichmentCsvView(View):
    def get(self, request, request_id):
        response = HttpResponse(content_type='text/csv')
        response["Content-Disposition"] = 'attachment; filename="analysis_enrichment.csv"'

        try:
            analysis_enrichment_csv(request_id, buffer=response)

            return response
        except AnalysisEnrichmentError:
            response.status_code = 404
            return response


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


class ListDownloadView(View):
    def head(self, request, list_name):
        files = chain(gene_lists_storage.listdir('.')[1], networks_storage.listdir('.')[1])

        name_regex = re.compile('^' + re.escape(list_name + os.path.extsep))

        if any(filter(name_regex.search, files)):
            return HttpResponse()
        else:
            return HttpResponseNotFound()

    def get(self, request, list_name):
        files = chain(gene_lists_storage.listdir('.')[1], networks_storage.listdir('.')[1])

        name_regex = re.compile('^' + re.escape(list_name + os.path.extsep))

        try:
            file = next(filter(name_regex.search, files))
            storage = next(filter(lambda s: os.path.isfile(s.path(file)), (gene_lists_storage, networks_storage)))

            if mimetypes.guess_type(file)[1] == 'gzip':
                return GzipFileResponse(storage.open(file, 'rb'), as_attachment=True)

            return FileResponse(storage.open(file, 'rb'), as_attachment=True)
        except (StopIteration, FileNotFoundError) as e:
            raise Http404('gene list not found') from e
