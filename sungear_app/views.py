import json
import logging
from collections import OrderedDict
from typing import Dict, List, Tuple

import pandas as pd
from django.core.cache import cache
from django.http import Http404, HttpResponseBadRequest, JsonResponse
from django.views import View
from sungear import SungearException, sungear

from querytgdb.models import Analysis
from querytgdb.utils import PandasJSONEncoder, clear_data, get_metadata

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

    df = pd.concat([df.iloc[:, ::2], df.iloc[:, 1::2]], axis=1)  # splitting odd and even rows

    if df.shape[1] < 2:
        raise SungearException("Sungear needs at least 2 analyses.")

    analyses = Analysis.objects.filter(
        pk__in=df.columns.get_level_values(1)
    )

    metadata = get_metadata(analyses)

    gene_lists = OrderedDict(
        (name, col.index[col.notna()]) for name, col in df.iteritems()
    )

    result, finished = sungear(gene_lists)

    return {
               **result,
               'metadata': [
                   (c, metadata.loc[c[1], :].to_dict()) for c in df.columns
               ]
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