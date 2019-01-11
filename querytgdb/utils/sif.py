import re
from collections import deque
from itertools import chain, groupby
from operator import methodcaller, itemgetter
from typing import Callable, Generator, Iterator, List, TextIO, Tuple, Union

import networkx as nx

WS = re.compile(r' +')


class InvalidEdge(ValueError):
    pass


Edge = Tuple[str, str, str]


def get_edges(s: Iterator[List[str]]) -> Generator[Edge, None, None]:
    """
    Build edge from a single line of a SIF file
    :param s:
    :return:
    """
    for line in s:
        if len(line) > 3:
            for e in line[:2]:
                yield (line[0], e, line[1])
        elif len(line) == 3:
            yield (line[0], line[2], line[1])
        else:
            raise InvalidEdge(f"A line needs to have either 1 or greater or equal than 3 columns: {line}")


def split_sif_lines(f: Iterator[str]) -> Iterator[List[str]]:
    """
    Splitting lines in a sif file

    :param f:
    :return:
    """
    text = deque()

    split_func: Callable[[str], List[str]]

    for line in f:
        text.append(line)
        if '\t' in line:
            split_func = methodcaller('split', '\t')
            break
    else:
        split_func = WS.split

    return map(split_func, map(methodcaller('rstrip', '\n'), filter(None, chain(text, f))))


def get_network(f: Union[TextIO, str]) -> nx.MultiDiGraph:
    """
    Get network from SIF file

    :param f:
    :return:
    """
    g = nx.MultiDiGraph()

    try:
        lines: List[str] = sorted(split_sif_lines(f.splitlines(keepends=True)), key=len)
    except AttributeError:
        lines = split_sif_lines(f)

    for key, group in groupby(lines, key=len):
        if key >= 3:
            g.add_edges_from(get_edges(group))
        elif key == 1:
            g.add_nodes_from(map(itemgetter(0), group))
        else:
            raise InvalidEdge(f"A line needs to have either 1 or greater or equal than 3 columns.")

    return g
