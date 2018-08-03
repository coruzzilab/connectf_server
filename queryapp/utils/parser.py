import re
import warnings
from collections import deque
from functools import partial, reduce
from operator import itemgetter
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union
from uuid import uuid4

import numpy as np
import pandas as pd
import pyparsing as pp
from django.db.models import Q
from django.db.utils import DatabaseError

from querytgdb.models import AnalysisIddata, Annotation, EdgeData, EdgeType, Interactions, MetaIddata, \
    ReferenceId, Regulation
from ..utils import cache_result, read_cached_result

__all__ = ['get_query_result', 'expand_ref_ids', 'ANNOTATIONS']


class DatabaseWarning(UserWarning):
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


def get_annotations() -> pd.DataFrame:
    try:
        anno = pd.DataFrame(
            Annotation.objects.values_list(
                'agi_id', 'ath_fullname', 'ath_gene_fam', 'ath_gene_type', 'ath_name', 'ath_id').iterator(),
            columns=['TARGET', 'Full Name', 'Gene Family', 'Type', 'Name', 'ath_id'])
        anno = anno.set_index('TARGET')
    except DatabaseError:
        warnings.warn(DatabaseWarning("No annotation data."))

        anno = pd.DataFrame(columns=['Full Name', 'Gene Family', 'Type', 'Name', 'ath_id'])
        anno.index.name = 'TARGET'

    return anno


ANNOTATIONS = get_annotations()

name = pp.Word(pp.alphanums + '-_.:')

point = pp.Literal('.')
e = pp.CaselessLiteral('e')
number = pp.Word(pp.nums)
plusorminus = pp.Literal('+') | pp.Literal('-')
integer = pp.Combine(pp.Optional(plusorminus) + number)

floatnum = pp.Combine(integer + pp.Optional(point + pp.Optional(number)) + pp.Optional(e + integer)).setParseAction(
    lambda toks: float(toks[0]))

and_ = pp.CaselessKeyword('and')
or_ = pp.CaselessKeyword('or')
andor = and_ | or_

opers = [
    (pp.CaselessKeyword("not"), 1, pp.opAssoc.RIGHT),
    (andor, 2, pp.opAssoc.LEFT)
]

quote_word = pp.Word(pp.alphanums + '-_. :')
single_quote = pp.Suppress(pp.Literal("'"))
double_quote = pp.Suppress(pp.Literal('"'))
quoted_name = double_quote + quote_word + double_quote | single_quote + quote_word + single_quote

modname = pp.Group((name | quoted_name)('key') + pp.oneOf('< = > >= <= !=')('oper') + (name | quoted_name)('value'))(
    'mod')

modifier = pp.Group(pp.Suppress('[') + pp.infixNotation(modname, opers) + pp.Suppress(']'))('modifier')

gene = name('gene_name')

expr = pp.infixNotation(gene, [(modifier, 1, pp.opAssoc.LEFT)] + opers)('query')


def is_name(key: str, item: Union[pp.ParseResults, Any]) -> bool:
    try:
        return item.getName() == key
    except AttributeError:
        return False


is_modifier = partial(is_name, 'modifier')
is_mod = partial(is_name, 'mod')

interactions = Interactions.objects.all()


def mod_to_str(curr: pp.ParseResults) -> str:
    if isinstance(curr, str):
        return curr
    else:
        return ' '.join(map(mod_to_str, curr))


def query_metadata(df: TargetFrame, key: str, value: str) -> pd.DataFrame:
    ref_ids = ReferenceId.objects.filter(
        Q(analysis_id_id__in=AnalysisIddata.objects.filter(
            analysis_type__iexact=key,
            analysis_value__iexact=value).values_list('analysis_id_id', flat=True)) |
        Q(meta_id_id__in=MetaIddata.objects.filter(
            meta_type__iexact=key,
            meta_value__iexact=value).values_list('meta_id_id', flat=True))).values_list('ref_id', flat=True)

    mask = pd.DataFrame(True, columns=df.columns, index=df.index)

    mask.loc[:, ~df.columns.get_level_values(1).isin(ref_ids)] = False

    return mask


def apply_comp_mod(df: TargetFrame, key: str, oper: str, value: Union[float, str]) -> pd.DataFrame:
    """
    apply Pvalue and Log2FC (fold change)
    """
    value = float(value)
    mask = pd.DataFrame(True, columns=df.columns, index=df.index)

    try:
        if oper == '=':
            return mask.where(df[(*df.name, key)] == value, False)
        elif oper == '>=':
            return mask.where(df[(*df.name, key)] >= value, False)
        elif oper == '<=':
            return mask.where(df[(*df.name, key)] <= value, False)
        elif oper == '>':
            return mask.where(df[(*df.name, key)] > value, False)
        elif oper == '<':
            return mask.where(df[(*df.name, key)] < value, False)
        elif oper == '!=':
            return mask.where(df[(*df.name, key)] != value, False)
        else:
            raise ValueError('invalid operator: {}'.format(oper))
    except KeyError:
        mask.loc[:, :] = False
        return mask


def apply_search_column(df: TargetFrame, key, value) -> pd.DataFrame:
    mask = pd.DataFrame(True, columns=df.columns, index=df.index)

    try:
        return mask.where(df[(*df.name, key)].str.contains(value, case=False, regex=False), False)
    except KeyError:
        mask.loc[:, :] = False
        return mask


COL_TRANSLATE = {
    'PVALUE': 'Pvalue',
    'FC': 'Log2FC'
}


def apply_has_column(df: TargetFrame, value) -> pd.DataFrame:
    try:
        value = COL_TRANSLATE[value]
    except KeyError:
        pass

    if (*df.name, value) in df:
        return pd.DataFrame(True, columns=df.columns, index=df.index)

    return pd.DataFrame(False, columns=df.columns, index=df.index)


def apply_has_edge_props(df: TargetFrame, value) -> pd.DataFrame:
    mask = pd.DataFrame(True, columns=df.columns, index=df.index)
    print(ReferenceId.objects.get(ref_id=df.name[1]))
    print(ANNOTATIONS.loc[df.index[df.notna().any(axis=1)], 'ath_id'])

    return mask


def get_mod(df: TargetFrame, query: Union[pp.ParseResults, pd.DataFrame]) -> pd.DataFrame:
    """
    Get ref_id from modifier to filter TF dataframe

    Careful not to modify original df
    """
    if isinstance(query, pp.ParseResults) and 'key' not in query:
        it = iter(query)
        stack = deque()

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
    elif isinstance(query, pp.ParseResults) and 'key' in query:
        key = query['key']
        oper = query['oper']
        value = query['value']

        if re.match(r'^pvalue$', key, flags=re.I):
            return df.groupby(level=[0, 1], axis=1).apply(apply_comp_mod, key='Pvalue', oper=oper, value=value)
        elif re.match(r'^fc$', key, flags=re.I):
            return df.groupby(level=[0, 1], axis=1).apply(apply_comp_mod, key='Log2FC', oper=oper, value=value)
        elif re.match(r'^edge$', key, flags=re.I):
            return df.groupby(level=[0, 1], axis=1).apply(apply_search_column, value=value, key='EDGE')
        # elif re.match(r'^edge_properties$', key, flags=re.I):
        #     return df.groupby(level=[0, 1], axis=1).apply(apply_has_edge_props, value=value)
        elif re.match(r'^has_column$', key, flags=re.I):
            value = value.upper()
            return df.groupby(level=[0, 1], axis=1).apply(apply_has_column, value=value)
        else:
            return query_metadata(df, key, value)
    else:
        return query


def get_tf_data(query: str, edges: Optional[List[str]] = None) -> TargetFrame:
    """
    Get data for single TF
    :param query:
    :param edges:
    :return:
    """
    df = TargetFrame(
        interactions.filter(db_tf_id__db_tf_agi__iexact=query).values_list(
            'edge_id__edge_name', 'target_id__agi_id', 'ref_id__ref_id').iterator(),
        columns=['EDGE', 'TARGET', 'REF'])
    if not df.empty:
        reg = TargetFrame(
            Regulation.objects.filter(ref_id__in=df['REF'].unique()).values_list(
                'ref_id', 'ath_id__agi_id', 'pvalue', 'foldchange').iterator(),
            columns=['REF', 'TARGET', 'Pvalue', 'Log2FC'])

        # somehow this is stored as strings for some reason
        reg['Pvalue'] = reg['Pvalue'].astype(float)
        reg['Log2FC'] = reg['Log2FC'].astype(float)

        df = df.merge(reg, on=['REF', 'TARGET'], how='left')

        if edges:
            anno = ANNOTATIONS['ath_id'].reset_index()

            edge_types = pd.DataFrame(EdgeType.objects.filter(name__in=edges).values_list('id', 'name').iterator(),
                                      columns=['edge_id', 'edge'])
            tf_ids = anno.loc[anno['TARGET'].str.contains(query, case=False, regex=False), 'ath_id']
            target_ids = anno.loc[anno['TARGET'].isin(df['TARGET'].unique()), 'ath_id']

            edge_data = pd.DataFrame(
                EdgeData.objects.filter(
                    tf_id__in=tf_ids,
                    target_id__in=target_ids
                ).values_list('tf_id', 'target_id', 'type_id').iterator(),
                columns=['source', 'target', 'edge_id']
            )

            edge_data = reduce(lambda x, y: x.str.cat(y, sep=',', na_rep='', join='outer'),
                               map(lambda x: x['edge'],
                                   map(itemgetter(1),
                                       (edge_data
                                        .merge(edge_types, on='edge_id')
                                        .drop('edge_id', axis=1)
                                        .set_index(['source', 'target'])
                                        .groupby('edge', squeeze=True))))
                               ).str.strip(',').reset_index()

            edge_data = edge_data.merge(anno, left_on='source', right_on='ath_id').merge(anno, left_on='target',
                                                                                         right_on='ath_id')
            edge_data = edge_data[['TARGET_x', 'TARGET_y', 'edge']]
            edge_data.columns = ['TF', 'TARGET', 'EDGE_PROPS']

            df = df.merge(edge_data.drop('TF', axis=1), on='TARGET')

        df = (df.pivot(index='TARGET', columns='REF')
              .swaplevel(0, 1, axis=1)
              .sort_index(axis=1, level=0, sort_remaining=False)
              .dropna(axis=1, how='all'))
    else:
        df = TargetFrame(columns=[(np.nan, 'EDGE')])

    df.columns = pd.MultiIndex.from_tuples((query, *c) for c in df.columns)
    df.filter_string += query

    return df


def get_all_tf(query: str, edges: Optional[List[str]] = None) -> TargetFrame:
    """
    Get data for all TFs at once
    :param query:
    :param edges:
    :return:
    """
    df = TargetFrame(
        interactions.values_list(
            'edge_id__edge_name', 'target_id__agi_id', 'ref_id__ref_id', 'db_tf_id__db_tf_agi').iterator(),
        columns=['EDGE', 'TARGET', 'REF', 'TF'])

    reg = TargetFrame(
        Regulation.objects.values_list(
            'ref_id', 'ath_id__agi_id', 'pvalue', 'foldchange').iterator(),
        columns=['REF', 'TARGET', 'Pvalue', 'Log2FC'])
    # somehow this is stored as strings for some reason
    reg['Pvalue'] = reg['Pvalue'].astype(float)
    reg['Log2FC'] = reg['Log2FC'].astype(float)

    df = df.merge(reg, on=['REF', 'TARGET'], how='left')

    if edges:
        anno = ANNOTATIONS['ath_id'].reset_index()
        edge_types = pd.DataFrame(EdgeType.objects.filter(name__in=edges).values_list('id', 'name').iterator(),
                                  columns=['edge_id', 'edge'])
        tf_ids = anno.loc[anno['TARGET'].isin(df['TF'].unique()), 'ath_id']
        target_ids = anno.loc[anno['TARGET'].isin(df['TARGET'].unique()), 'ath_id']

        edge_data = pd.DataFrame(
            EdgeData.objects.filter(
                tf_id__in=tf_ids,
                target_id__in=target_ids
            ).values_list('tf_id', 'target_id', 'type_id').iterator(),
            columns=['source', 'target', 'edge_id']
        )

        edge_data = reduce(lambda x, y: x.str.cat(y, sep=',', na_rep='', join='outer'),
                           map(lambda x: x['edge'],
                               map(itemgetter(1),
                                   (edge_data
                                    .merge(edge_types, on='edge_id')
                                    .drop('edge_id', axis=1)
                                    .set_index(['source', 'target'])
                                    .groupby('edge', squeeze=True))))
                           ).str.strip(',').reset_index()

        edge_data = edge_data.merge(anno, left_on='source', right_on='ath_id').merge(anno, left_on='target',
                                                                                     right_on='ath_id')
        edge_data = edge_data[['TARGET_x', 'TARGET_y', 'edge']]
        edge_data.columns = ['TF', 'TARGET', 'EDGE_PROPS']

        df = df.merge(edge_data, on=['TF', 'TARGET'], how='left')

    df = (df.set_index(['TF', 'REF', 'TARGET'])
          .unstack(level=[0, 1])
          .reorder_levels([1, 2, 0], axis=1)
          .sort_index(axis=1, level=[0, 1], sort_remaining=False)
          .dropna(how='all', axis=1))

    if query == 'oralltfs':
        df.filter_string += 'oralltfs'
    elif query == 'andalltfs':
        df = df[df.loc[:, (slice(None), slice(None), 'EDGE')].notna().all(axis=1)]
        df.filter_string += 'andalltfs'
    else:
        raise ValueError('invalid query')

    return df


def get_suffix(prec: TargetFrame, succ: TargetFrame) -> Tuple[str, str]:
    if prec.filter_string == succ.filter_string:
        return '_' + str(uuid4()), '_' + str(uuid4())

    return ' "' + prec.filter_string + '"', ' "' + succ.filter_string + '"'


def get_tf(query: Union[pp.ParseResults, str, TargetFrame], edges: Optional[List[str]] = None) -> TargetFrame:
    """
    Query TF DataFrame according to query
    :param query:
    :param edges:
    :return:
    """
    if isinstance(query, pp.ParseResults):
        it = iter(query)
        stack = deque()

        try:
            while True:
                curr = next(it)
                if curr in ('and', 'or'):
                    prec, succ = get_tf(stack.pop(), edges), get_tf(next(it), edges)

                    filter_string = prec.filter_string
                    if curr == 'and':
                        filter_string += ' and '

                        if prec.include and succ.include:
                            df = prec.merge(succ, how='inner', left_index=True, right_index=True,
                                            suffixes=get_suffix(prec, succ))
                        elif not prec.include and succ.include:
                            df = succ.loc[~succ.index.isin(prec.index), :]
                        elif prec.include and not succ.include:
                            df = prec.loc[~prec.index.isin(succ.index), :]
                        else:  # not prec.include and not succ.include
                            df = prec.merge(succ, how='outer', left_index=True, right_index=True,
                                            suffixes=get_suffix(prec, succ))
                            df.include = False
                    else:
                        filter_string += ' or '

                        # doesn't make much sense using not with or, but oh well
                        if prec.include and succ.include:
                            df = prec.merge(succ, how='outer', left_index=True, right_index=True,
                                            suffixes=get_suffix(prec, succ))
                        elif not prec.include and succ.include:
                            df = succ
                        elif prec.include and not succ.include:
                            df = prec
                        else:
                            df = prec.merge(succ, how='inner', left_index=True, right_index=True,
                                            suffixes=get_suffix(prec, succ))
                            df.include = False
                    filter_string += succ.filter_string
                    df = df.groupby(level=[0, 1], axis=1).filter(lambda x: x.notna().any(axis=None))
                    df.filter_string = filter_string
                    stack.append(df)
                elif curr == 'not':
                    succ = get_tf(next(it), edges)
                    succ.include = not succ.include
                    succ.filter_string = 'not ' + succ.filter_string
                    stack.append(succ)
                elif is_modifier(curr):
                    prec = get_tf(stack.pop(), edges)
                    mod = get_mod(prec, curr)
                    prec = prec[mod].dropna(how='all')

                    # filter out empty tfs
                    prec = prec.groupby(level=[0, 1], axis=1).filter(lambda x: x.notna().any(axis=None))
                    prec.filter_string += '[' + mod_to_str(curr) + ']'

                    stack.append(prec)
                else:
                    stack.append(curr)
        except StopIteration:
            return get_tf(stack.pop(), edges)
    elif isinstance(query, (TargetFrame, TargetSeries)):
        return query
    else:
        if query in {'andalltfs', 'oralltfs'}:
            return get_all_tf(query, edges)

        return get_tf_data(query, edges)


def reorder_data(df: TargetFrame) -> TargetFrame:
    """
    Order by TF with most edges, then analysis with most edges within tf
    :param df:
    :return:
    """
    analysis_order = df.loc[:, (slice(None), slice(None), slice(None), 'EDGE')].count(
        axis=1, level=2).sum().sort_values(ascending=False)
    meta_order = df.loc[:, (slice(None), slice(None), slice(None), 'EDGE')].count(axis=1, level=1).sum().sort_values(
        ascending=False)
    tf_order = df.loc[:, (slice(None), slice(None), slice(None), 'EDGE')].count(axis=1, level=0).sum().sort_values(
        ascending=False)

    df = (df.reindex(labels=analysis_order.index, axis=1, level=2)
          .reindex(labels=meta_order.index, axis=1, level=1)
          .reindex(labels=tf_order.index, axis=1, level=0))
    return df


def get_metadata(ids: Sequence) -> TargetFrame:
    refs = ReferenceId.objects.filter(ref_id__in=ids)
    df = pd.concat([
        TargetFrame(refs.values_list(
            'ref_id',
            'analysis_id__analysisiddata__analysis_type',
            'analysis_id__analysisiddata__analysis_value').iterator(),
                    columns=['REF', 'KEY', 'VALUE']),
        TargetFrame(
            refs.values_list(
                'ref_id',
                'meta_id__metaiddata__meta_type',
                'meta_id__metaiddata__meta_value').iterator(),
            columns=['REF', 'KEY', 'VALUE'])
    ], ignore_index=True)

    df = df.dropna(how='all', subset=['KEY', 'VALUE'])
    df = df.set_index(['REF', 'KEY'])
    df = df.unstack(level=0)
    df.columns = df.columns.droplevel(level=0)

    return df


def get_tf_count(df: TargetFrame) -> TargetSeries:
    counts = df.loc[:, (slice(None), slice(None), slice(None), 'EDGE')].count(axis=1)
    counts.name = 'TF Count'

    return counts


def expand_ref_ids(df: TargetFrame, level: Optional[Union[str, int]] = None) -> TargetFrame:
    df = df.copy()

    if level is None:
        ref_ids = df.columns
    else:
        ref_ids = df.columns.levels[level]

    full_ids = {ref_id: (meta_id, analysis_id) for ref_id, meta_id, analysis_id in
                ReferenceId.objects.filter(
                    ref_id__in=ref_ids).values_list('ref_id',
                                                    'meta_id__meta_fullid',
                                                    'analysis_id__analysis_fullid').iterator()}

    if level is None:
        df.columns = pd.MultiIndex.from_tuples(full_ids[ref_id] for ref_id in df.columns)
    elif isinstance(df.columns, pd.MultiIndex):
        df.columns = pd.MultiIndex.from_tuples((*c[:level], *full_ids[c[level]], *c[level + 1:]) for c in df.columns)
    else:
        raise ValueError('Please specify level to expand.')

    return df


def parse_query(query: str, edges: Optional[List[str]] = None) -> TargetFrame:
    parse = expr.parseString(query, parseAll=True)

    result = get_tf(parse.get('query'), edges)

    if result.empty or not result.include:
        raise ValueError('empty query')

    return result


def trim_edges(df: TargetFrame) -> TargetFrame:
    df.loc[:, (slice(None), slice(None), slice(None), 'EDGE')] = \
        df.loc[:, (slice(None), slice(None), slice(None), 'EDGE')].apply(
            lambda x: x.str.replace('.+(?=INDUCED|REPRESSED|REGULATED|BOUND)', '',
                                    regex=True, flags=re.I))

    return df


def induce_repress_count(df: TargetFrame):
    return TargetSeries((df[(*df.name, 'EDGE')].str.contains('induced', case=False).sum(),
                         df[(*df.name, 'EDGE')].str.contains('repressed', case=False).sum()),
                        index=['induced', 'repressed'])


def get_stats(result: pd.DataFrame) -> Dict[str, Any]:
    return {
        'total': result.loc[:, (slice(None), slice(None), slice(None), 'EDGE')].groupby(
            level=[0, 1, 2], axis=1).count().sum()
    }


def get_query_result(query: Optional[str] = None,
                     user_lists: Optional[Tuple[pd.DataFrame, Dict]] = None,
                     edges: Optional[List[str]] = None,
                     cache_path: Optional[Union[str, Path]] = None) -> Tuple[pd.DataFrame, pd.DataFrame, Dict]:
    if query is None and cache_path is None:
        raise ValueError("Need query or cache_path")

    if query is not None:
        result = parse_query(query, edges)
        metadata = get_metadata(result.columns.get_level_values(1))
        result = expand_ref_ids(result, 1)

        stats = get_stats(result)

        if cache_path:
            result.to_pickle(cache_path + '/tabular_output_unfiltered.pickle.gz')

        if user_lists:
            result = result[result.index.isin(user_lists[0].index)]

        if cache_path:  # cache here
            result.to_pickle(cache_path + '/tabular_output.pickle.gz')
            metadata.to_pickle(cache_path + '/metadata.pickle.gz')

            if user_lists:
                cache_result(user_lists, cache_path + '/target_genes.pickle.gz')
    else:
        result = pd.read_pickle(cache_path + '/tabular_output.pickle.gz')
        metadata = pd.read_pickle(cache_path + '/metadata.pickle.gz')
        stats = get_stats(pd.read_pickle(cache_path + '/tabular_output_unfiltered.pickle.gz'))

        try:
            user_lists = read_cached_result(cache_path + '/target_genes.pickle.gz')
        except FileNotFoundError:
            pass

    result = reorder_data(result)

    result = trim_edges(result)

    counts = get_tf_count(result)

    result = pd.concat([counts, result], axis=1)

    if user_lists:
        result = user_lists[0].merge(result, left_index=True, right_index=True, how='inner')
        result = result.sort_values(['User List Count', 'User List'])
    else:
        result = pd.concat([
            TargetFrame(np.nan, columns=['User List', 'User List Count'], index=result.index),
            result
        ], axis=1)

    result = ANNOTATIONS.drop('ath_id', axis=1).merge(result, how='right', left_index=True, right_index=True)

    result = result.sort_values('TF Count', ascending=False, kind='mergesort')

    # return statistics here as well
    return result, metadata, stats
