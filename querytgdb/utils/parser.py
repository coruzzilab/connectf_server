import logging
import operator
import re
from collections import UserDict, defaultdict, deque
from concurrent.futures import ThreadPoolExecutor
from functools import partial, reduce
from operator import and_, itemgetter, methodcaller, or_
from typing import Any, Callable, Deque, Dict, List, Optional, Sequence, Tuple, Union
from uuid import UUID, uuid4

import numpy as np
import pandas as pd
import pyparsing as pp
from django.conf import settings
from django.core.cache import cache, caches
from django.core.exceptions import MultipleObjectsReturned, ObjectDoesNotExist
from django.db.models import Q

from querytgdb.models import Analysis, Annotation, EdgeData, EdgeType, Interaction, Regulation
from querytgdb.utils import async_loader
from ..utils import CaselessDict, clear_data, get_metadata as get_meta_df
from ..utils.file import UserGeneLists

logger = logging.getLogger(__name__)

mem_cache = caches['mem']

NAMED_QUERIES = CaselessDict(getattr(settings, 'NAMED_QUERIES', {}))

PVALUE = 'Pvalue'
LOG2FC = 'Log2FC'


class QueryError(ValueError):
    pass


class TargetFrame(pd.DataFrame):
    _metadata = ['include', 'filter_string']

    @property
    def _constructor(self):
        return TargetFrame

    @property
    def _constructor_sliced(self):
        return TargetSeries

    def __init__(self, *args, include=True, **kwargs):
        super().__init__(*args, **kwargs)
        self.include = include
        self.filter_string = ''


class TargetSeries(pd.Series):
    @property
    def _constructor(self):
        return TargetSeries

    @property
    def _constructor_expanddim(self):
        return TargetFrame


name = pp.Word(pp.pyparsing_unicode.alphanums + '-_.:/\\')
pl_name = pp.CaselessKeyword('$filter_tf')('filter_tf')  # placeholder for expansion

and_oper = pp.CaselessKeyword('and')
or_oper = pp.CaselessKeyword('or')

opers = [
    (pp.CaselessKeyword("not"), 1, pp.opAssoc.RIGHT),
    (and_oper | or_oper, 2, pp.opAssoc.LEFT)
]

mod_comp = pp.oneOf('< = > >= <= != <>')

quoted_name = pp.QuotedString('"', escChar='\\') | pp.QuotedString("'", escChar='\\')

pvalue = pp.CaselessKeyword('pvalue')
log2fc = pp.CaselessKeyword('log2fc')
id_ = pp.CaselessKeyword('id')

number_keywords = (pvalue | log2fc | id_).setParseAction(lambda toks: toks[0].lower())

additional_edge = pp.CaselessKeyword('additional_edge')
targeted_by = pp.CaselessKeyword('targeted_by')

string_keywords = (additional_edge | targeted_by).setParseAction(lambda toks: toks[0].lower())

modname = pp.Group(
    (number_keywords('key') + mod_comp('oper') + pp.pyparsing_common.number('value')) |
    (string_keywords('key') + mod_comp('oper') + (name | quoted_name)('value')) |
    ((name | quoted_name)('key') + mod_comp('oper') + (name | quoted_name)('value'))
)('mod')

modifier = pp.Group(pp.Suppress('[') + pp.infixNotation(modname, opers) + pp.Suppress(']'))('modifier')

column_filter = pp.Group(pp.Suppress('{') + pp.delimitedList(modname, ',') + pp.Suppress('}'))('column_filter')

expr = pp.Forward()

all_tfs = pp.CaselessKeyword('all_tfs')
multitype = pp.CaselessKeyword('multitype')


def parse_prebuilt_query(tocs):
    return expr.parseString(NAMED_QUERIES[tocs[0]], parseAll=True)


named_query = reduce(lambda a, b: a | b, map(pp.CaselessKeyword, NAMED_QUERIES.keys())).setParseAction(
    parse_prebuilt_query)

gene = (all_tfs | multitype | named_query | name)('gene_name')

# additional functions

func_name = pp.Word(pp.pyparsing_unicode.alphas.lower(), pp.pyparsing_unicode.alphanums)('func_name')
query_func = (func_name + pp.Suppress('(') + pp.delimitedList(quoted_name | pp.Word(pp.alphanums + '$_ '),
                                                              ',') + pp.Suppress(')'))('function')

expr <<= pp.infixNotation(gene, [(modifier | column_filter, 1, pp.opAssoc.LEFT)] + opers)('query') ^ query_func


def is_name(key: str, item: Union[pp.ParseResults, Any]) -> bool:
    try:
        return item.getName() == key
    except AttributeError:
        return False


is_modifier = partial(is_name, 'modifier')
is_mod = partial(is_name, 'mod')
is_column_filter = partial(is_name, 'column_filter')


def mod_to_str(curr: pp.ParseResults) -> str:
    if isinstance(curr, str):
        return curr

    if isinstance(curr, pp.ParseResults) and curr.getName() != 'mod':
        return '(' + ' '.join(map(mod_to_str, curr)) + ')'

    try:
        return ' '.join(map(mod_to_str, curr))
    except TypeError:
        return str(curr)


def query_metadata(df: TargetFrame, key: str, value: str) -> pd.DataFrame:
    ref_ids = Analysis.objects.filter(
        (Q(analysisdata__key__name__iexact=key) & Q(analysisdata__value__iexact=value))
    ).values_list('id', flat=True)

    if not ref_ids:
        return pd.DataFrame(False, columns=df.columns, index=df.index)

    mask = pd.DataFrame(True, columns=df.columns, index=df.index)

    mask.loc[:, ~df.columns.get_level_values(1).isin(ref_ids)] = False

    return mask


OPERS = {
    '>': operator.gt,
    '>=': operator.ge,
    '<': operator.lt,
    '<=': operator.le,
    '=': operator.eq,
    '!=': operator.ne,
    '<>': operator.ne
}


def apply_comp_mod(df: TargetFrame, key: str, oper: str, value: float) -> pd.DataFrame:
    """
    apply Pvalue and Log2FC (fold change)
    """
    try:
        c = df[(*df.name, key)]
    except KeyError:
        return pd.DataFrame(False, columns=df.columns, index=df.index)

    mask = pd.DataFrame(True, columns=df.columns, index=df.index)

    try:
        return mask.where(OPERS[oper](c, value), False)
    except KeyError as e:
        raise ValueError('invalid operator: {}'.format(oper)) from e


def apply_search_column(df: TargetFrame, key, value) -> pd.DataFrame:
    mask = pd.DataFrame(True, columns=df.columns, index=df.index)

    try:
        return mask.where(df[(*df.name, key)].str.contains(value, case=False, regex=False), False)
    except KeyError:
        mask.loc[:, :] = False
        return mask


def match_id(df: TargetFrame, oper: str, analysis_id: Union[str, int]):
    """
    Filter dataframe by analysis_id
    :param df:
    :param oper:
    :param analysis_id:
    :return:
    """
    mask = pd.DataFrame(False, columns=df.columns, index=df.index)
    mask.loc[:, (slice(None), OPERS[oper](df.columns.get_level_values(1), int(analysis_id)), slice(None))] = True

    return mask


COL_TRANSLATE = {
    'PVALUE': PVALUE,
    'LOG2FC': LOG2FC,
    'ADDITIONAL_EDGE': 'ADD_EDGES'
}


def apply_has_column(df: TargetFrame, value) -> pd.DataFrame:
    try:
        value = COL_TRANSLATE[value]
    except KeyError:
        pass

    if (*df.name, value) in df:
        return pd.DataFrame(True, columns=df.columns, index=df.index)

    return pd.DataFrame(False, columns=df.columns, index=df.index)


def apply_has_add_edges(df: TargetFrame, analyses, anno_ids, value) -> pd.DataFrame:
    mask = pd.DataFrame(True, columns=df.columns, index=df.index)
    try:
        edge_type = EdgeType.objects.get(name__iexact=value)
        target_ids = set(EdgeData.objects.filter(
            type=edge_type,
            tf_id=analyses.get(pk=df.name[1]).tf_id
        ).values_list('target_id', flat=True).iterator())

        target_ids &= set(anno_ids.loc[df.index[df.notna().any(axis=1)]].unique())

        mask.loc[~df.index.isin(anno_ids.index[anno_ids.isin(target_ids)]), :] = False
    except (ObjectDoesNotExist, MultipleObjectsReturned):
        mask.loc[:, :] = False

    return mask


def match_targeted_by(df, oper, value):
    frac = False
    if '%' in value:
        value = float(value.rstrip('% ')) / 100
        frac = True
    else:
        value = float(value)

    if 0 <= value < 1:
        frac = True

    mask = pd.DataFrame(False, columns=df.columns, index=df.index)

    cleared = df.pipe(clear_data)

    targeted_count = cleared.count(axis=1)

    op_func = OPERS[oper]

    if frac:
        rows = op_func(targeted_count / cleared.shape[1], value)
    else:
        rows = op_func(targeted_count, value)

    mask.loc[rows, :] = True

    return mask


def get_mod(df: TargetFrame, query: Union[pp.ParseResults, pd.DataFrame]) -> pd.DataFrame:
    """
    Get ref_id from modifier to filter TF dataframe

    Careful not to modify original df
    """
    if isinstance(query, pp.ParseResults):
        if 'key' not in query:
            it = iter(query)
            stack: Deque[Union[pp.ParseResults, pd.DataFrame]] = deque()

            try:
                while True:
                    curr = next(it)
                    if curr in ('and', 'or'):
                        prec, succ = get_mod(df, stack.pop()), get_mod(df, next(it))
                        if curr == 'and':
                            stack.append(prec & succ)
                        else:
                            stack.append(prec | succ)
                    elif curr == 'not':
                        succ = get_mod(df, next(it))
                        stack.append(~succ)
                    else:
                        stack.append(curr)
            except StopIteration:
                return get_mod(df, stack.pop())
        else:
            key, oper, value = itemgetter('key', 'oper', 'value')(query)
            if key == 'pvalue':
                return df.groupby(level=[0, 1], axis=1).apply(apply_comp_mod, key=PVALUE, oper=oper, value=value)
            elif key == 'log2fc':
                return df.groupby(level=[0, 1], axis=1).apply(apply_comp_mod, key=LOG2FC, oper=oper, value=value)
            elif key == 'additional_edge':
                analyses = Analysis.objects.filter(pk__in=df.columns.get_level_values(1)).prefetch_related('tf')
                anno_ids = async_loader['annotations'].loc[df.index, 'id']
                return df.groupby(level=[0, 1], axis=1).apply(apply_has_add_edges,
                                                              analyses=analyses,
                                                              anno_ids=anno_ids,
                                                              value=value)
            elif key == 'id':
                return match_id(df, oper, value)
            elif key == 'targeted_by':
                return match_targeted_by(df, oper, value)
            else:
                return query_metadata(df, key, value)
    return query


def gene_to_ids(metadata, ids):
    d: Dict[str, set] = defaultdict(set)
    for idx in ids:
        d[metadata.at[idx, 'gene_id']].add(idx)

    return d


def get_column_filter(df: TargetFrame, filter_list) -> TargetFrame:
    metadata = get_meta_df(Analysis.objects.filter(pk__in=df.columns.get_level_values(1)))
    metadata.columns = metadata.columns.str.lower()
    result = []
    for query in filter_list:
        key, oper, value = itemgetter('key', 'oper', 'value')(query)
        key = key.lower()

        if key == 'pvalue':
            c = OPERS[oper](df.loc[:, (slice(None), slice(None), [PVALUE])], value)
            result.append(c.columns[c.any(axis=0)].get_level_values(1))
        elif key == 'log2fc':
            c = OPERS[oper](df.loc[:, (slice(None), slice(None), [LOG2FC])], value)
            result.append(c.columns[c.any(axis=0)].get_level_values(1))
        elif key in metadata.columns:
            result.append(metadata.index[metadata[key] == value])
        else:
            raise QueryError(f"Invalid column filter: {key}")

    result = list(map(partial(gene_to_ids, metadata), result))
    valid_genes = reduce(and_, map(methodcaller('keys'), result))
    if not valid_genes:
        return TargetFrame(columns=pd.MultiIndex(levels=[[], [], []]))

    return df.loc[:, df.columns.get_level_values(1).isin(reduce(or_, (r[g] for r in result for g in valid_genes)))]


def add_edges(df: pd.DataFrame, edges: List[str]) -> pd.DataFrame:
    """
    Add additional edge data of query to result dataframe
    :param df:
    :param edges:
    :return:
    """
    anno = async_loader['annotations']['id'].reset_index()

    edge_types = pd.DataFrame(
        EdgeType.objects.filter(name__in=edges).values_list('id', 'name', 'directional').iterator(),
        columns=['edge_id', 'edge', 'directional'])

    try:
        tf_ids = anno.loc[anno['TARGET'].isin(df['TF'].unique()), 'id']
    except KeyError:
        tf_ids = Annotation.objects.filter(analysis__in=df['ANALYSIS'].unique()).values_list('pk', flat=True)

    edge_data = pd.DataFrame(
        EdgeData.objects.filter(
            tf_id__in=tf_ids
        ).values_list('tf_id', 'target_id', 'type_id').iterator(),
        columns=['source', 'target', 'edge_id']
    )

    edge_data = edge_data.loc[edge_data['target'].isin(df['id']), :]

    edge_data = (edge_data
                 .merge(edge_types[['edge_id', 'edge']], on='edge_id')
                 .drop('edge_id', axis=1)
                 .set_index(['source', 'target']))

    if edge_data.empty:
        raise ValueError("No Edge Data")

    edge_data = pd.concat(map(itemgetter(1), edge_data.groupby('edge')), axis=1)

    row_num, col_num = edge_data.shape

    if col_num > 1:
        edge_data = edge_data.iloc[:, 0].str.cat(
            map(itemgetter(1), edge_data.iloc[:, 1:].iteritems()),
            sep=',', na_rep='', join='inner').str.strip(',')
    else:
        edge_data = edge_data.fillna('')

    edge_data = edge_data.reset_index()

    edge_data = edge_data.merge(anno, left_on='source', right_on='id')
    edge_data = edge_data[['TARGET', 'target', 'edge']]
    edge_data.columns = ['TF', 'id', 'ADD_EDGES']

    if 'TF' in df:
        return df.merge(edge_data, on=['TF', 'id'], how='left')
    return df.merge(edge_data.drop('TF', axis=1), on='id', how='left')


def initialize_column_name(col_name) -> Tuple[str, str, str]:
    return col_name, "", str(uuid4())


def get_tf_data(query: str,
                edges: Optional[List[str]] = None,
                tf_filter_list: Optional[pd.Series] = None,
                target_filter_list: Optional[pd.Series] = None) -> TargetFrame:
    """
    Get data for single TF
    :param query:
    :param edges:
    :param tf_filter_list:
    :param target_filter_list:
    :return:
    """
    anno = async_loader['annotations']

    if (tf_filter_list is not None and tf_filter_list.str.contains(rf'^{re.escape(query)}$', flags=re.I).any()) \
            or tf_filter_list is None:
        analyses = Analysis.objects.filter(tf__gene_id__iexact=query)

        if not analyses.exists():
            raise ValueError(f'"{query}" is not in database')

        df = TargetFrame(
            Interaction.objects.filter(analysis__in=analyses).values_list(
                'target_id', 'analysis_id').iterator(),
            columns=['id', 'ANALYSIS'])
        df = df.merge(anno['id'].reset_index(), on=['id'])
        if target_filter_list is not None:
            df = df[df['TARGET'].str.upper().isin(target_filter_list.str.upper())]
        df = df.reindex(columns=['TARGET', 'ANALYSIS', 'id'])
    else:
        analyses = []
        df = TargetFrame(columns=['TARGET', 'ANALYSIS', 'id'])

    if not df.empty:
        reg = TargetFrame(
            Regulation.objects.filter(analysis__in=analyses).values_list(
                'analysis_id', 'target_id', 'p_value', 'foldchange').iterator(),
            columns=['ANALYSIS', 'id', PVALUE, LOG2FC])

        df.insert(2, 'EDGE', '+')

        expressions = Analysis.objects.filter(
            pk__in=analyses,
            analysisdata__key__name='EXPERIMENT_TYPE',
            analysisdata__value__iexact='expression'
        ).values_list('pk', flat=True)

        df.loc[df['ANALYSIS'].isin(expressions), 'EDGE'] = '*'

        if not reg.empty:
            df = df.merge(reg, on=['ANALYSIS', 'id'], how='left')
            df.loc[~(df[LOG2FC].isna() & df[PVALUE].isna()), 'EDGE'] = np.nan

        if edges:
            try:
                df = add_edges(df, edges)
            except ValueError:
                pass

        df = df.drop('id', axis=1)

        df = (df.pivot(index='TARGET', columns='ANALYSIS')
              .swaplevel(0, 1, axis=1)
              .sort_index(axis=1, level=0, sort_remaining=False)
              .dropna(axis=1, how='all'))
    else:
        df = TargetFrame(columns=[(np.nan, 'EDGE')])

    try:
        query = anno.index[np.argmax(anno['id'] == analyses[0].tf_id)].upper()
    except IndexError:
        query = query.upper()

    q = initialize_column_name(query)
    df.columns = pd.MultiIndex.from_tuples((q, *c) for c in df.columns)
    df.filter_string = query

    return df


def get_all_interaction(qs):
    return TargetFrame(qs.iterator(), columns=['id', 'ANALYSIS'])


def get_all_regulation(qs):
    return TargetFrame(qs.iterator(), columns=['ANALYSIS', 'id', PVALUE, LOG2FC])


def get_all_analyses(qs):
    return TargetFrame(qs.iterator(), columns=['ANALYSIS', 'TF'])


def get_all_df(query: str,
               tf_filter_list: Optional[pd.Series] = None,
               target_filter_list: Optional[pd.Series] = None) -> TargetFrame:
    qs = Interaction.objects.values_list('target_id', 'analysis_id')
    anno = async_loader['annotations']

    if tf_filter_list is not None:
        qs = qs.filter(analysis__tf_id__in=anno.loc[anno.index.str.upper().isin(tf_filter_list.str.upper()), 'id'])

    with ThreadPoolExecutor(max_workers=3) as executor:
        interaction_task = executor.submit(get_all_interaction, qs)
        regulation_task = executor.submit(get_all_regulation,
                                          Regulation.objects.values_list('analysis_id', 'target_id', 'p_value',
                                                                         'foldchange'))
        analysis_task = executor.submit(get_all_analyses, Analysis.objects.values_list('id', 'tf__gene_id'))

        df = interaction_task.result()

        df = df.merge(anno['id'].reset_index(), on='id')
        df = df.reindex(columns=['TARGET', 'ANALYSIS', 'id'])

        if query == "multitype":
            a = pd.DataFrame(Analysis.objects.filter(
                analysisdata__key__name__iexact="EXPERIMENT_TYPE"
            ).values_list('id', 'tf_id', 'analysisdata__value', named=True).iterator())

            a = a.groupby('tf_id').filter(lambda x: x['analysisdata__value'].nunique() > 1)

            df = df[df['ANALYSIS'].isin(a['id'])]

        if target_filter_list is not None:
            df = df[df['TARGET'].str.upper().isin(target_filter_list.str.upper())]

        analyses = analysis_task.result()
        df = df.merge(analyses, on='ANALYSIS')

        reg = regulation_task.result()

    df = df.merge(reg, on=['ANALYSIS', 'id'], how='left')

    # additional restrictions here as well
    if query == 'andalltfs':
        df = df[df.loc[:, (slice(None), slice(None), ['EDGE', 'Log2FC'])].notna().all(axis=1)]

    return df


def replace_filter_str(col: Tuple[str, str, str], filter_string: str) -> Tuple[str, str, str]:
    return col[0], filter_string, col[2]


def get_all_tf(query: str,
               edges: Optional[List[str]] = None,
               tf_filter_list: Optional[pd.Series] = None,
               target_filter_list: Optional[pd.Series] = None) -> TargetFrame:
    """
    Get data for all TFs at once
    :param query:
    :param edges:
    :param tf_filter_list:
    :param target_filter_list:
    :return:
    """
    if tf_filter_list is None and target_filter_list is None:
        df = mem_cache.get_or_set(query, partial(get_all_df, query))
    else:
        df = get_all_df(query, tf_filter_list, target_filter_list)

    if df.empty:
        raise ValueError("No data in database.")

    if edges:
        try:
            df = add_edges(df, edges)
        except ValueError:
            pass

    df = df.drop('id', axis=1)

    df.insert(3, 'EDGE', np.nan)
    df['EDGE'] = df['EDGE'].where(df[LOG2FC].notna() | df[PVALUE].notna(), '+')

    expressions = Analysis.objects.filter(
        analysisdata__key__name='EXPERIMENT_TYPE',
        analysisdata__value__iexact='expression'
    ).values_list('pk', flat=True)

    df.loc[df[LOG2FC].isna() & df[PVALUE].isna() & (df['ANALYSIS'].isin(expressions)), 'EDGE'] = '*'

    df = (df.set_index(['TF', 'ANALYSIS', 'TARGET'])
          .unstack(level=[0, 1])
          .reorder_levels([1, 2, 0], axis=1)
          .sort_index(axis=1, level=[0, 1], sort_remaining=False)
          .dropna(how='all', axis=1))

    col_names = {c: initialize_column_name(c) for c in df.columns.levels[0]}

    df = df.rename(columns=col_names, level=0)

    df.filter_string += query

    return df


def get_tf(query: Union[pp.ParseResults, str, TargetFrame],
           edges: Optional[List[str]] = None,
           tf_filter_list: Optional[pd.Series] = None,
           target_filter_list: Optional[pd.Series] = None) -> TargetFrame:
    """
    Query TF DataFrame according to query
    :param query:
    :param edges:
    :param tf_filter_list:
    :param target_filter_list:
    :return:
    """
    if isinstance(query, pp.ParseResults):
        it = iter(query)
        stack: Deque[Union[pd.DataFrame, str, pp.ParseResults]] = deque()

        try:
            while True:
                curr = next(it)
                if curr in ('and', 'or'):
                    prec, succ = get_tf(stack.pop(), edges, tf_filter_list, target_filter_list), \
                                 get_tf(next(it), edges, tf_filter_list, target_filter_list)

                    filter_string = prec.filter_string
                    if curr == 'and':
                        filter_string += ' and '

                        if prec.include and succ.include:
                            df = prec.merge(succ, how='inner', left_index=True, right_index=True)
                        elif not prec.include and succ.include:
                            df = succ.loc[~succ.index.isin(prec.index), :]
                        elif prec.include and not succ.include:
                            df = prec.loc[~prec.index.isin(succ.index), :]
                        else:  # not prec.include and not succ.include
                            df = prec.merge(succ, how='outer', left_index=True, right_index=True)
                            df.include = False
                    else:
                        filter_string += ' or '

                        # doesn't make much sense using not with or, but oh well
                        if prec.include and succ.include:
                            df = prec.merge(succ, how='outer', left_index=True, right_index=True)
                        elif not prec.include and succ.include:
                            df = succ
                        elif prec.include and not succ.include:
                            df = prec
                        else:
                            df = prec.merge(succ, how='inner', left_index=True, right_index=True)
                            df.include = False
                    filter_string = '(' + filter_string + succ.filter_string + ')'

                    try:
                        df = df.dropna(axis=1, how='all')
                    except IndexError:
                        # beware of the shape of indices and columns
                        df = TargetFrame(columns=pd.MultiIndex(levels=[[], [], []]))

                    df.filter_string = filter_string

                    if curr == 'and':
                        # df = df.rename(columns=partial(replace_filter_str, filter_string=filter_string), level=0)
                        df.columns = pd.MultiIndex.from_tuples([(replace_filter_str(c[0], filter_string), c[1], c[2]) for c in df.columns])

                    stack.append(df)

                elif curr == 'not':
                    succ = get_tf(next(it), edges, tf_filter_list, target_filter_list)
                    succ.include = not succ.include
                    succ.filter_string = 'not ' + succ.filter_string
                    stack.append(succ)
                elif is_modifier(curr):
                    prec = get_tf(stack.pop(), edges, tf_filter_list, target_filter_list)
                    mod = get_mod(prec, curr)
                    prec = prec[mod].dropna(how='all').dropna(how='all', axis=1)  # filter out empty tfs

                    prec.filter_string += f'[{mod_to_str(curr[0])}]'
                    prec = prec.rename(columns=partial(replace_filter_str, filter_string=prec.filter_string), level=0)

                    stack.append(prec)
                elif is_column_filter(curr):
                    prec = get_tf(stack.pop(), edges, tf_filter_list, target_filter_list)
                    prec = get_column_filter(prec, curr)

                    stack.append(prec)
                else:
                    stack.append(curr)
        except StopIteration:
            return get_tf(stack.pop(), edges, tf_filter_list, target_filter_list)
    elif isinstance(query, (TargetFrame, TargetSeries)):
        return query
    elif isinstance(query, str):
        if query.lower() in {'andalltfs', 'all_tfs', 'multitype'}:
            return get_all_tf(query.lower(), edges, tf_filter_list, target_filter_list)

        return get_tf_data(query, edges, tf_filter_list, target_filter_list)
    else:
        raise ValueError(query)


def reorder_data(df: TargetFrame) -> TargetFrame:
    """
    Order by TF with most edges, then analysis with most edges within tf
    :param df:
    :return:
    """
    edges = clear_data(df)

    analysis_order = edges.groupby(level=1, axis=1).count().sum().sort_values(ascending=False)
    tf_order = edges.groupby(level=0, axis=1).count().sum()
    tf_total = tf_order.groupby(by=list(map(itemgetter(0), tf_order.index))).sum()
    tf_reorder = sorted(tf_order.index, key=lambda i: (tf_total[i[0]], tf_order.at[i]), reverse=True)

    return (df.reindex(columns=analysis_order.index, level=1)
            .reindex(columns=tf_reorder, level=0))


def get_metadata(ids: Sequence) -> pd.DataFrame:
    df = get_meta_df(ids)
    df = df.T
    df.columns = df.columns.rename('ANALYSIS')
    df.index = df.index.rename('KEY')

    return df


def get_tf_count(df: TargetFrame) -> TargetSeries:
    counts = df.loc[:, (slice(None), slice(None), ['EDGE', 'Log2FC'])].count(axis=1)
    counts.name = 'Edge Count'

    return counts


def add_tf_count(df: TargetFrame) -> TargetFrame:
    counts = df.pipe(get_tf_count)

    df = df.copy()
    df.columns = df.columns.to_flat_index()
    df.insert(0, 'Edge Count', counts)

    return df


def expand_ref_ids(df: pd.DataFrame, level: Optional[Union[str, int]] = None) -> pd.DataFrame:
    df = df.copy()

    if level is None:
        analysis_ids = df.columns
    else:
        analysis_ids = df.columns.levels[level]

    full_ids = {a: a.name for a in
                Analysis.objects.filter(
                    pk__in=analysis_ids
                ).prefetch_related('analysisdata_set', 'analysisdata_set__key')}

    if level is None:
        df = df.rename(columns=full_ids)
    elif isinstance(df.columns, pd.MultiIndex):
        df = df.rename(columns=full_ids, level=level)
    else:
        raise ValueError('Please specify level to expand.')

    return df


class QueryFuncNotFound(KeyError, QueryError):
    pass


class QueryFuncs(UserDict):
    def __getitem__(self, item):
        try:
            return super().__getitem__(item)
        except KeyError as e:
            raise QueryFuncNotFound(e) from e

    def register(self, func: Callable):
        self[func.__name__] = func
        return func


QUERY_FUNCS = QueryFuncs()


@QUERY_FUNCS.register
def expand(query: str, oper: str, *args, tf_filter_list: Optional[pd.Series] = None, **kwargs):
    if tf_filter_list is None or tf_filter_list.empty:
        raise QueryError('Filter TF list required')

    result = f" {oper} ".join(re.sub(r'\$filter_tf', t, f"({query})", flags=re.I) for t in tf_filter_list)

    return parse_query(result, tf_filter_list=tf_filter_list, **kwargs)


def parse_query(query: str,
                edges: Optional[List[str]] = None,
                tf_filter_list: Optional[pd.Series] = None,
                target_filter_list: Optional[pd.Series] = None) -> TargetFrame:
    try:
        parse = expr.parseString(query, parseAll=True)

        if parse.getName() == 'function':
            fname, *args = parse
            return QUERY_FUNCS[fname](*args,
                                      edges=edges,
                                      tf_filter_list=tf_filter_list,
                                      target_filter_list=target_filter_list)

        result = get_tf(parse.get('query'), edges, tf_filter_list, target_filter_list)

        if result.empty or not result.include:
            raise QueryError('empty query')

        result.columns = result.columns.set_levels(result.columns.levels[1].astype(int), level=1)  # ensure integer type

        result = reorder_data(result)

        return result
    except pp.ParseException as e:
        raise QueryError("Could not parse query") from e


def get_total(df: pd.DataFrame) -> pd.DataFrame:
    return df.pipe(clear_data).groupby(level=[0, 1], axis=1).count().sum()


def induce_repress_count(result: TargetFrame) -> pd.DataFrame:
    fc_result = result.loc[:, (slice(None), slice(None), ['Log2FC'])]
    return pd.DataFrame([(fc_result > 0).sum(axis=0), (fc_result < 0).sum(axis=0)], index=['induced', 'repressed'])


Id = Tuple[Tuple[str, str, str], int]
Properties = Dict[str, Any]
Ids = Dict[Id, Properties]


def get_result_ids(df: pd.DataFrame) -> Ids:
    return {c: {'name': f'{c[0][0]}_{c[0][1]}_{c[1]}', 'show': True, 'version': 0} for c in df.columns.droplevel(2)}


def filter_df_by_ids(df: pd.DataFrame, ids: Ids) -> pd.DataFrame:
    show_ids = {k for k, v in ids.items() if v['show']}

    if not show_ids:
        raise ValueError("All analyses hidden")

    if len(show_ids) == len(ids):
        return df

    return df.loc[:, list(map(lambda c: itemgetter(0, 1)(c) in show_ids, df.columns))]


def get_query_result(query: Optional[str] = None,
                     uid: Optional[Union[str, UUID]] = None,
                     user_lists: Optional[UserGeneLists] = None,
                     tf_filter_list: Optional[pd.Series] = None,
                     target_filter_list: Optional[pd.Series] = None,
                     edges: Optional[List[str]] = None,
                     size_limit: Optional[int] = None) \
        -> Tuple[pd.DataFrame, pd.DataFrame, Dict, Union[str, UUID], Ids]:
    """
    Get query result from query string or cache

    :param query:
    :param uid:
    :param user_lists:
    :param tf_filter_list:
    :param target_filter_list:
    :param edges:
    :param size_limit:
    :return:
    """
    if uid is None:
        uid = uuid4()

    if query is not None:  # Check if cached
        result = parse_query(query, edges, tf_filter_list, target_filter_list)
        metadata = get_metadata(result.columns.get_level_values(1))
        ids = get_result_ids(result)

        cache.set_many({
            f'{uid}/tabular_output_unfiltered': result,
            f'{uid}/metadata': metadata,
            f'{uid}/analysis_ids': ids
        })
    else:
        data = cache.get_many([
            f'{uid}/tabular_output_unfiltered',
            f'{uid}/analysis_ids'
        ])

        result, ids = itemgetter(
            f'{uid}/tabular_output_unfiltered',
            f'{uid}/analysis_ids'
        )(data)

        result = filter_df_by_ids(result, ids)
        metadata = get_metadata(result.columns.get_level_values(1))

        cache.set(f'{uid}/metadata', metadata)

    stats = {
        'total': get_total(result)
    }

    if user_lists is None and query is None:
        user_lists = cache.get(f'{uid}/target_genes')

    if user_lists is not None:
        result = result[result.index.str.upper().isin(user_lists[0].index.str.upper())].dropna(axis=1, how='all')

        if result.empty:
            raise QueryError("Empty result (user list too restrictive).")

        result = reorder_data(result)  # reorder again here due to filtering

        stats['edge_counts'] = get_total(result)
    else:
        stats['edge_counts'] = stats['total']

    stats['induce_repress_count'] = result.pipe(induce_repress_count)

    cache.set_many({f'{uid}/tabular_output': result})

    logger.info(f"Unfiltered Dataframe size: {result.size}")

    if size_limit is not None and result.size > size_limit:
        raise QueryError("Result too large.")

    result = result.pipe(add_tf_count)

    if user_lists:
        result = user_lists[0].merge(result,
                                     left_on=user_lists[0].index.str.upper(),
                                     right_on=result.index.str.upper(), how='inner')
        result = result.rename(columns={'key_0': 'TARGET'}).set_index('TARGET')
        result = result.sort_values(['User List Count', 'User List'])
    else:
        result.insert(0, 'User List Count', np.nan)
        result.insert(0, 'User List', np.nan)

    result = async_loader['annotations'].drop('id', axis=1).merge(result, how='right', left_index=True,
                                                                  right_index=True)

    result = result.sort_values('Edge Count', ascending=False, kind='mergesort')

    logger.info(f"Dataframe size: {result.size}")

    # return statistics here as well
    return result, metadata, stats, uid, ids
