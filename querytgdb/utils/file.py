import gzip
import mimetypes
import os.path
import re
from collections import OrderedDict
from contextlib import closing
from io import TextIOWrapper
from typing import Dict, Hashable, IO, Optional, Set, TextIO, Tuple

import numpy as np
import pandas as pd
from django.core.exceptions import SuspiciousFileOperation
from django.core.files.storage import Storage
from django.http.request import HttpRequest
from pandas.errors import EmptyDataError, ParserError


class BadFile(ValueError):
    pass


class BadNetwork(BadFile):
    pass


def get_file(request: HttpRequest, key: Hashable, storage: Optional[Storage] = None) -> Tuple[
    Optional[TextIO], Optional[str]]:
    """
    Get file or file name from the request.

    :param request:
    :param key:
    :param storage:
    :return:
    """
    if request.FILES and key in request.FILES:
        return TextIOWrapper(request.FILES[key]), 'post'
    elif key in request.POST and storage is not None:
        try:
            name = request.POST[key]
            name_regex = re.compile('^' + re.escape(name + os.path.extsep))

            directories, files = storage.listdir('.')

            file = next(filter(name_regex.search, files))

            if mimetypes.guess_type(file)[1] == 'gzip':
                return gzip.open(storage.path(file), 'rt'), 'storage'

            return storage.open(file, 'r'), 'storage'
        except (FileNotFoundError, SuspiciousFileOperation, StopIteration):
            pass

    return None, None


def gene_list_to_df(gene_to_name: Dict[str, Set[str]]) -> pd.DataFrame:
    return pd.DataFrame(
        ((key, ', '.join(val), len(val)) for key, val in gene_to_name.items()),
        columns=['TARGET', 'User List', 'User List Count']
    ).set_index('TARGET')


def get_gene_lists(f: TextIO) -> Tuple[pd.DataFrame, Dict[str, Set[str]]]:
    """
    Get gene lists from the uploaded target genes file.

    :param f:
    :return:
    """
    gene_to_name: Dict[str, Set[str]] = OrderedDict()
    name_to_gene: Dict[str, Set[str]] = OrderedDict()

    with closing(f) as gene_file:
        list_name = 'default_list'
        for line in gene_file:
            line = line.strip()
            if line.startswith('>'):
                list_name = line.lstrip('>').strip()
            else:
                gene_to_name.setdefault(line, set()).add(list_name)
                name_to_gene.setdefault(list_name, set()).add(line.upper())

    if not gene_to_name:
        raise BadFile("Target Gene list empty")

    df = gene_list_to_df(gene_to_name)

    return df, name_to_gene


def get_genes(f: TextIO) -> pd.Series:
    with closing(f) as g:
        s = pd.Series(g.readlines())

    if s.empty:
        raise BadFile("Filter TF list is empty")

    s = s.str.strip()
    s = s[~(s.str.startswith('>') | s.str.startswith(';'))].reset_index(drop=True)

    return s


Network = Tuple[str, pd.DataFrame]

NETWORK_MSG = "Network must have source, edge, target columns. Can have an additional forth column of scores."


def get_network(f: IO) -> Network:
    """
    Parse uploaded file into dataframe
    :param f:
    :return:
    """
    try:
        df = pd.read_csv(f, delim_whitespace=True, header='infer')
    except (ParserError, UnicodeDecodeError, EmptyDataError) as e:
        raise BadNetwork(NETWORK_MSG) from e

    name = os.path.basename(getattr(f, 'name', 'default'))

    rows, cols = df.shape

    if cols == 2:
        df.columns = ['source', 'target']
    elif cols == 3:
        if np.issubdtype(df.dtypes[2], np.number):  # use last column as score if number
            df.columns = ['source', 'target', 'score']
        else:
            df.columns = ['source', 'edge', 'target']
    elif cols == 4:
        df.columns = ['source', 'edge', 'target', 'score']
    else:
        raise BadNetwork(NETWORK_MSG)

    if 'score' in df:
        df['rank'] = df['score'].rank(method='max', ascending=False)
        df = df.sort_values('rank')
    else:
        df['rank'] = np.arange(1, df.shape[0] + 1)

    if 'edge' not in df:
        df.insert(1, 'edge', name)

    return name, df


def network_to_lists(network: Network) -> Tuple[pd.DataFrame, Dict[str, Set[str]]]:
    """
    Makes network into user_lists format
    :param network:
    :return:
    """
    name, data = network

    gene_to_name: Dict[str, Set[str]] = OrderedDict()
    name_to_gene: Dict[str, Set[str]] = OrderedDict()

    genes = data[['source', 'target']].stack().unique()

    for g in genes:
        gene_to_name.setdefault(g, set()).add(name)

    name_to_gene[name] = set(genes)

    df = gene_list_to_df(gene_to_name)

    return df, name_to_gene


def merge_network_lists(user_lists: Tuple[pd.DataFrame, OrderedDict],
                        network: Tuple[str, pd.DataFrame]) -> Tuple[pd.DataFrame, Dict[str, Set[str]]]:
    """
    Make graphs into user_lists format and merging with user_lists
    :param user_lists:
    :param network:
    :return:
    """
    graph_lists = network_to_lists(network)

    name_to_gene = graph_lists[1]

    for k, v in user_lists[1].items():
        name_to_gene.setdefault(k, set()).update(v)

    df = graph_lists[0].merge(user_lists[0], left_index=True, right_index=True, how='outer')
    names = df["User List_x"].str.cat(df["User List_y"], sep=', ', na_rep='').str.strip(', ').rename("User List")
    count = (df["User List Count_x"].fillna(0) + df["User List Count_y"].fillna(0)).astype(int).rename(
        "User List Count")

    df = pd.concat([names, count], axis=1)

    return df, name_to_gene


def network_to_filter_tfs(network: Tuple[str, pd.DataFrame]) -> pd.Series:
    """
    Use source nodes (and isolated nodes) as filter tfs for dataframe

    :param network:
    :return:
    """
    return pd.Series(network[1]['source'].unique())
