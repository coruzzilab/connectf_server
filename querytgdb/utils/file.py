import re
from collections import OrderedDict, defaultdict
from contextlib import closing
from io import TextIOWrapper
from itertools import chain, groupby
from operator import itemgetter, methodcaller
from typing import DefaultDict, Dict, Generator, Hashable, Iterable, Iterator, List, Optional, Set, TextIO, Tuple, Union

import networkx as nx
import pandas as pd
from django.core.exceptions import SuspiciousFileOperation
from django.core.files.storage import Storage
from django.http.request import HttpRequest

Graphs = DefaultDict[str, nx.DiGraph]

WS = re.compile(r' +')


def get_file(request: HttpRequest, key: Hashable, storage: Optional[Storage] = None) -> Optional[TextIO]:
    """
    Get file or file name from the request.

    :param request:
    :param key:
    :param storage:
    :return:
    """
    if request.FILES and key in request.FILES:
        return TextIOWrapper(request.FILES[key])
    elif key in request.POST and storage is not None:
        try:
            return storage.open("{}.txt".format(request.POST[key]), 'r')
        except (FileNotFoundError, SuspiciousFileOperation):
            pass


def gene_list_to_df(gene_to_name: Dict[str, Set[str]]) -> pd.DataFrame:
    return pd.DataFrame(
        ((key, ', '.join(val), len(val)) for key, val in gene_to_name.items()),
        columns=['TARGET', 'User List', 'User List Count']
    ).set_index('TARGET')


def get_gene_lists(f: TextIO) -> Tuple[pd.DataFrame, OrderedDict]:
    """
    Get gene lists from the uploaded target genes file.

    :param f:
    :return:
    """
    gene_to_name = OrderedDict()
    name_to_gene = OrderedDict()

    with closing(f) as gene_file:
        list_name = 'default_list'
        for line in gene_file:
            line = line.strip()
            if line.startswith('>'):
                list_name = line.lstrip('>').strip()
            else:
                gene_to_name.setdefault(line, set()).add(list_name)
                name_to_gene.setdefault(list_name, set()).add(line)

    df = gene_list_to_df(gene_to_name)

    return df, name_to_gene


def get_genes(f: TextIO) -> pd.Series:
    with closing(f) as g:
        s = pd.Series(g.readlines())
    s = s.str.strip()
    s = s[~(s.str.startswith('>') | s.str.startswith(';'))].reset_index(drop=True)

    return s


def split_sif_lines(f: TextIO) -> Iterator[List[str]]:
    with closing(f) as g:
        text = g.read()

    if '\t' in text:
        func = methodcaller('split', '\t')
    else:
        func = WS.split

    return map(func, filter(None, text.splitlines()))


Edge = Tuple[str, str, str, int]


def build_edge(lines: Iterable[List[str]]) -> Generator[Union[Edge, str], None, None]:
    """
    Yields source, target, edge, rank tuple

    Discards all lines with out edges

    :param lines:
    :return:
    """
    i = 1
    for line in lines:
        if len(line) > 3:
            for e in line[:2]:
                yield (line[0], e, line[1], i)
                i += 1
        elif len(line) == 3:
            yield (line[0], line[2], line[1], i)
            i += 1


def build_node(lines: List[List[str]]) -> Iterable[str]:
    """
    Return node name from sif file for nodes with no edges

    :param lines:
    :return:
    """
    return map(itemgetter(0), filter(lambda l: len(l) < 3, lines))


def get_network(f: TextIO) -> Graphs:
    """
    Parse uploaded file into Graphs
    :param f:
    :return:
    """
    graphs = defaultdict(nx.DiGraph)

    lines = list(split_sif_lines(f))

    edges = groupby(sorted(build_edge(filter(lambda l: len(l) > 2, lines)), key=itemgetter(2, 3), reverse=True),
                    key=itemgetter(2))
    nodes = list(build_node(lines))

    for name, group in edges:
        graphs[name].add_edges_from((e[0], e[1], {'name': name, 'rank': e[3]}) for e in group)
        graphs[name].add_nodes_from(nodes)

    if not graphs:
        graphs[f"default_graph"].add_nodes_from(nodes)

    return graphs


def network_to_lists(graphs: Graphs) -> Tuple[pd.DataFrame, OrderedDict]:
    """
    Makes graphs into user_lists format
    :param graphs:
    :return:
    """
    gene_to_name = OrderedDict()
    name_to_gene = OrderedDict()

    for name, graph in graphs.items():
        name_to_gene[name] = set(graph.nodes)
        for n in graph.nodes:
            gene_to_name.setdefault(n, set()).add(name)

    df = gene_list_to_df(gene_to_name)

    return df, name_to_gene


def merge_network_lists(user_lists: Tuple[pd.DataFrame, OrderedDict],
                        graphs: Graphs) -> Tuple[pd.DataFrame, OrderedDict]:
    """
    Make graphs into user_lists format and merging with user_lists
    :param user_lists:
    :param graphs:
    :return:
    """
    graph_lists = network_to_lists(graphs)

    name_to_gene = graph_lists[1]

    for k, v in user_lists[1].items():
        name_to_gene.setdefault(k, set()).update(v)

    df = graph_lists[0].merge(user_lists[0], left_index=True, right_index=True, how='outer')
    names = df["User List_x"].str.cat(df["User List_y"], sep=', ', na_rep='').str.strip(', ').rename("User List")
    count = (df["User List Count_x"].fillna(0) + df["User List Count_y"].fillna(0)).astype(int).rename(
        "User List Count")

    df = pd.concat([names, count], axis=1)

    return df, name_to_gene


def network_to_filter_tfs(graphs: Graphs) -> pd.Series:
    """
    Use source nodes (and isolated nodes) as filter tfs for dataframe

    :param graphs:
    :return:
    """
    return pd.Series(list(
        chain(
            chain.from_iterable((n for n, o in g.out_degree if o > 0) for g in graphs.values()),
            chain.from_iterable((n for n, d in g.degree if d == 0) for g in graphs.values())
        )
    ))
