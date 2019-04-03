import os.path
from collections import OrderedDict, defaultdict
from functools import reduce
from itertools import groupby
from operator import itemgetter, or_
from typing import DefaultDict, Dict, List
from uuid import uuid4

import numpy as np
import pandas as pd
from django.conf import settings
from django.core.files.storage import FileSystemStorage
from django.http import Http404, HttpResponse
from django.views import View

from querytgdb.models import Analysis
from querytgdb.utils import clear_data

static_storage = FileSystemStorage(settings.QUERY_CACHE)


class SungearException(ValueError):
    pass


def get_vertices(n: int, r: float = 0.4375, center=(0.5, 0.5), start=0) -> np.ndarray:
    return r * np.array([
        [np.cos(i - start), np.sin(i - start)]
        for i in np.linspace(2 * np.pi, 0, n, endpoint=False)
    ]) + center


# Create your views here.
class SungearView(View):
    def get(self, request, request_id):
        try:
            cache_path = static_storage.path(f"{request_id}_pickle")

            df = pd.read_pickle(os.path.join(cache_path, 'tabular_output.pickle.gz'))
            df = clear_data(df)

            if df.shape[1] < 2:
                raise SungearException("Sungear needs at least 2 analyses.")

            analyses = Analysis.objects.filter(
                pk__in=df.columns.get_level_values(1)
            ).prefetch_related('analysisdata_set', 'tf')

            gene_lists = OrderedDict(
                (analyses.get(pk=name[1]).name, set(col.index[col.notna()])) for name, col in df.iteritems()
            )

            total_genes = reduce(or_, gene_lists.values())
            gene_to_list: DefaultDict[str, List] = defaultdict(list)

            for name, l in gene_lists.items():
                for g in (total_genes & l):
                    gene_to_list[g].append(name)

            num_lists = len(gene_lists)

            vertices = pd.DataFrame(get_vertices(num_lists))
            vertices.index = gene_lists.keys()

            print(vertices)

            # if num_lists > 15:
            #     raise SungearException("Too many analyses (>15)")

            intersects: Dict[str, List] = OrderedDict()

            for name, group in groupby(
                    sorted(((tuple(sorted(v)), k) for k, v in gene_to_list.items()),
                           key=lambda x: (len(x[0]), x[0], x[1])),
                    key=itemgetter(0)):
                uid = str(uuid4())
                if len(name) > 1:
                    intersects[uid] = [name, vertices.loc[list(name), :].mean(axis=0), tuple(map(itemgetter(1), group))]

            print(intersects)

            return HttpResponse(request_id)
        except FileNotFoundError as e:
            raise Http404 from e
