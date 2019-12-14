import json
import logging
from collections import OrderedDict
from functools import reduce
from operator import or_
from typing import Dict, List, Tuple

import pandas as pd
from django.core.cache import cache
from django.http import Http404, HttpResponseBadRequest, JsonResponse
from django.views import View
from sungear import SungearException, sungear

from querytgdb.models import Analysis
from querytgdb.utils import PandasJSONEncoder, async_loader, clear_data, get_metadata

logger = logging.getLogger(__name__)


class SungearNotFound(SungearException):
    pass


class FilterListNotFound(SungearException):
    pass


def get_sungear(uid, filter_genes: List[str] = None) -> Tuple[Dict, bool]:
    df = cache.get(f'{uid}/tabular_output')
    if df is None:
        raise SungearNotFound()

    df = clear_data(df)

    if filter_genes:
        df = df.loc[df.index.isin(filter_genes), :].dropna(how='all', axis=1)
    if df.shape[1] < 2:
        raise SungearException("Sungear needs at least 2 analyses.")

    df = pd.concat([df.iloc[:, ::2], df.iloc[:, 1::2]], axis=1)  # splitting odd and even rows

    analyses = Analysis.objects.filter(
        pk__in=df.columns.get_level_values(1)
    )

    metadata = get_metadata(analyses)

    gene_lists = OrderedDict(
        (name, col.index[col.notna()]) for name, col in df.iteritems()
    )

    genes = list(reduce(or_, gene_lists.values()))
    anno = async_loader['annotations'].loc[genes, ['Full Name', 'Name']].rename(
        columns={'Full Name': 'name', 'Name': 'symbol'})
    anno = anno[(anno['name'] != "") | (anno['symbol'] != "")]

    result, finished = sungear(gene_lists)

    return {
               **result,
               'metadata': [
                   (c, metadata.loc[c[1], :].to_dict()) for c in df.columns
               ],
               'gene_annotation': anno.to_dict('index')
           }, finished and not filter_genes


# Create your views here.
class SungearView(View):
    def get(self, request, request_id):
        try:
            res = cache.get(f'{request_id}/sungear')
            if res is None:
                res, can_cache = get_sungear(request_id)
                if can_cache:
                    cache.set(f'{request_id}/sungear', res)
            return JsonResponse(res, encoder=PandasJSONEncoder)
        except SungearNotFound:
            raise Http404
        except SungearException as e:
            return HttpResponseBadRequest(e)

    def post(self, request, request_id):
        try:
            try:
                genes = json.loads(request.body)['genes']

                if not genes:
                    raise FilterListNotFound('empty genes')

                res, can_cache = get_sungear(request_id, genes)
            except (json.JSONDecodeError, KeyError, FilterListNotFound):
                res = cache.get(f'{request_id}/sungear')
                if res is None:
                    res, can_cache = get_sungear(request_id)
                    if can_cache:
                        cache.set(f'{request_id}/sungear', res)
            return JsonResponse(res, encoder=PandasJSONEncoder)
        except SungearNotFound:
            raise Http404
        except SungearException as e:
            return HttpResponseBadRequest(e)
