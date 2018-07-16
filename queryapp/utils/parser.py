import re
from collections import deque
from functools import partial
from pathlib import Path
from typing import Any, Dict, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
import pyparsing as pp
from django.db.models import Q
from pandas.api.extensions import register_dataframe_accessor

from querytgdb.models import AnalysisIddata, Annotation, DAPdata, Interactions, MetaIddata, ReferenceId, Regulation

__all__ = ['get_query_result', 'expand_ref_ids']


@register_dataframe_accessor('target')
class TargetAccessor:
    def __init__(self, pandas_obj):
        self._obj = pandas_obj
        self._include = True

    @property
    def include(self):
        return self._include

    @include.setter
    def include(self, value):
        self._include = value


def get_annotations() -> pd.DataFrame:
    anno = pd.DataFrame(
        Annotation.objects.values_list(
            'agi_id', 'ath_fullname', 'ath_gene_fam', 'ath_gene_type', 'ath_name').iterator(),
        columns=['TARGET', 'Full Name', 'Gene Family', 'Type', 'Name'])
    anno = anno.set_index('TARGET')

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


def query_metadata(df: pd.DataFrame, key: str, value: str):
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


def apply_comp_mod(df: pd.DataFrame, key: str, oper: str, value: Union[float, str]) -> pd.DataFrame:
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


def apply_search_column(df: pd.DataFrame, key, value) -> pd.DataFrame:
    mask = pd.DataFrame(True, columns=df.columns, index=df.index)

    return mask.where(df[(*df.name, key)].str.contains(value, case=False, regex=False), False)


def get_mod(df: pd.DataFrame, query: pp.ParseResults):
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
        elif re.match(r'^dap$', key, flags=re.I):
            return df.groupby(level=[0, 1], axis=1).apply(apply_search_column, value=value, key='DAP')
        else:
            return query_metadata(df, key, value)
    else:
        return query


def get_tf_data(query: str) -> pd.DataFrame:
    """
    Get data for single TF
    :param query:
    :return:
    """
    df = pd.DataFrame(
        interactions.filter(db_tf_id__db_tf_agi__exact=query).values_list(
            'edge_id__edge_name', 'target_id__agi_id', 'ref_id__ref_id').iterator(),
        columns=['EDGE', 'TARGET', 'REF'])
    if df.empty:
        raise ValueError('No Data. (invalid query tf)')

    reg = pd.DataFrame(
        Regulation.objects.filter(ref_id__in=df['REF'].unique()).values_list(
            'ref_id', 'ath_id__agi_id', 'pvalue', 'foldchange').iterator(),
        columns=['REF', 'TARGET', 'Pvalue', 'Log2FC'])

    # somehow this is stored as strings for some reason
    reg['Pvalue'] = reg['Pvalue'].astype(float)
    reg['Log2FC'] = reg['Log2FC'].astype(float)

    dap = pd.Series(
        DAPdata.objects.filter(db_tfid__db_tf_agi__exact=query).values_list('ath_id__agi_id', flat=True).iterator())
    df = df.merge(reg, on=['REF', 'TARGET'], how='left')
    if not dap.empty:
        df.loc[df['TARGET'].isin(dap), 'DAP'] = 'Present'
    df = df.pivot(index='TARGET', columns='REF')
    df = df.swaplevel(0, 1, axis=1)
    df = df.sort_index(axis=1, level=0, sort_remaining=False)

    df = df.dropna(axis=1, how='all')

    df.columns = pd.MultiIndex.from_tuples((query, *c) for c in df.columns)

    return df


def get_all_tf(query: str) -> pd.DataFrame:
    """
    Get data for all TFs at once
    :param query:
    :return:
    """
    df = pd.DataFrame(
        interactions.values_list(
            'edge_id__edge_name', 'target_id__agi_id', 'ref_id__ref_id', 'db_tf_id__db_tf_agi').iterator(),
        columns=['EDGE', 'TARGET', 'REF', 'TF'])

    reg = pd.DataFrame(
        Regulation.objects.values_list(
            'ref_id', 'ath_id__agi_id', 'pvalue', 'foldchange').iterator(),
        columns=['REF', 'TARGET', 'Pvalue', 'Log2FC'])
    # somehow this is stored as strings for some reason
    reg['Pvalue'] = reg['Pvalue'].astype(float)
    reg['Log2FC'] = reg['Log2FC'].astype(float)

    dap = pd.DataFrame(
        DAPdata.objects.values_list('db_tfid__db_tf_agi', 'ath_id__agi_id').iterator(),
        columns=['TF', 'TARGET'])
    dap['DAP'] = 'Present'
    df = df.merge(reg, on=['REF', 'TARGET'], how='left')
    df = df.merge(dap, on=['TF', 'TARGET'], how='left')
    df = df.set_index(['TF', 'REF', 'TARGET']).unstack(level=[0, 1])
    df = df.reorder_levels([1, 2, 0], axis=1)
    df = df.sort_index(axis=1, level=[0, 1], sort_remaining=False)
    df = df.dropna(how='all', axis=1)

    if query == 'oralltf':
        pass
    elif query == 'andalltf':
        df = df[df.loc[:, (slice(None), slice(None), 'EDGE')].notna().all(axis=1)]
    else:
        raise ValueError('invalid query')

    return df


def get_tf(query: Union[pp.ParseResults, str, pd.DataFrame]) -> pd.DataFrame:
    """
    Query TF DataFrame according to query
    :param query:
    :return:
    """
    if isinstance(query, pp.ParseResults):
        it = iter(query)
        stack = deque()

        try:
            while True:
                curr = next(it)
                if curr in ('and', 'or'):
                    prec, succ = get_tf(stack.pop()), get_tf(next(it))

                    if curr == 'and':
                        if prec.target.include and succ.target.include:
                            df = prec.merge(succ, how='inner', left_index=True, right_index=True)
                        elif not prec.target.include and succ.target.include:
                            df = succ.loc[~succ.index.isin(prec.index), :]
                        elif prec.target.include and not succ.target.include:
                            df = prec.loc[~prec.index.isin(succ.index), :]
                        else:  # not prec.target.include and not succ.target.include
                            df = prec.merge(succ, how='outer', left_index=True, right_index=True)
                            df.target.include = False
                    else:
                        # doesn't make much sense using not with or, but oh well
                        if prec.target.include and succ.target.include:
                            df = prec.merge(succ, how='outer', left_index=True, right_index=True)
                        elif not prec.target.include and succ.target.include:
                            df = succ
                        elif prec.target.include and not succ.target.include:
                            df = prec
                        else:
                            df = prec.merge(succ, how='inner', left_index=True, right_index=True)
                            df.target.include = False
                    inc = df.target.include
                    df = df.groupby(level=[0, 1], axis=1).filter(lambda x: x.notna().any(axis=None))
                    df.target.include = inc
                    stack.append(df)
                elif curr == 'not':
                    succ = get_tf(next(it))
                    succ.target.include = not succ.target.include
                    stack.append(succ)
                elif is_modifier(curr):
                    prec = get_tf(stack.pop())
                    mod = get_mod(prec, curr)
                    prec = prec[mod].dropna(how='all')

                    # filter out empty tfs
                    prec = prec.groupby(level=[0, 1], axis=1).filter(lambda x: x.notna().any(axis=None))

                    stack.append(prec)
                else:
                    stack.append(curr)
        except StopIteration:
            return get_tf(stack.pop())
    elif isinstance(query, (pd.DataFrame, pd.Series)):
        return query
    else:
        if query in ('andalltf', 'oralltf'):
            return get_all_tf(query)

        return get_tf_data(query)


def reorder_data(df: pd.DataFrame) -> pd.DataFrame:
    """
    Order by TF with most edges, then analysis with most edges within tf
    :param df:
    :return:
    """
    # expand the reference section into analysis and meta names here
    analysis_order = df.loc[:, (slice(None), slice(None), slice(None), 'EDGE')].count(
        axis=1, level=2).sum().sort_values(ascending=False)
    meta_order = df.loc[:, (slice(None), slice(None), slice(None), 'EDGE')].count(axis=1, level=1).sum().sort_values(
        ascending=False)
    tf_order = df.loc[:, (slice(None), slice(None), slice(None), 'EDGE')].count(axis=1, level=0).sum().sort_values(
        ascending=False)

    df = df.reindex(labels=analysis_order.index, axis=1, level=2)
    df = df.reindex(labels=meta_order.index, axis=1, level=1)
    df = df.reindex(labels=tf_order.index, axis=1, level=0)
    return df


def get_metadata(ids: Sequence) -> pd.DataFrame:
    refs = ReferenceId.objects.filter(ref_id__in=ids)
    df = pd.DataFrame(refs.values_list(
        'ref_id',
        'analysis_id__analysisiddata__analysis_type',
        'analysis_id__analysisiddata__analysis_value').iterator(),
                      columns=['REF', 'KEY', 'VALUE'])

    df = pd.concat([
        df,
        pd.DataFrame(
            refs.values_list(
                'ref_id',
                'meta_id__metaiddata__meta_type',
                'meta_id__metaiddata__meta_value').iterator(),
            columns=['REF', 'KEY', 'VALUE'])
    ])

    df = df.set_index(['REF', 'KEY'])
    df = df.unstack(level=0)
    df.columns = df.columns.droplevel(level=0)

    return df


def get_tf_count(df: pd.DataFrame) -> pd.Series:
    counts = df.loc[:, (slice(None), slice(None), slice(None), 'EDGE')].count(axis=1)
    counts.name = 'TF Count'

    return counts


def expand_ref_ids(df: pd.DataFrame, level: Optional[Union[str, int]] = None) -> pd.DataFrame:
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


def parse_query(query: str) -> pd.DataFrame:
    parse = expr.parseString(query, parseAll=True)

    result = get_tf(parse.get('query'))

    if result.empty or not result.target.include:
        raise ValueError('empty query')

    return result


def trim_edges(df: pd.DataFrame) -> pd.DataFrame:
    df.loc[:, (slice(None), slice(None), slice(None), 'EDGE')] = \
        df.loc[:, (slice(None), slice(None), slice(None), 'EDGE')].apply(
            lambda x: x.str.replace('.+(?=INDUCED|REPRESSED|REGULATED|BOUND)', '',
                                    regex=True, flags=re.I))

    return df


def induce_repress_count(df: pd.DataFrame):
    return pd.Series((df[(*df.name, 'EDGE')].str.contains('induced', case=False).sum(),
                      df[(*df.name, 'EDGE')].str.contains('repressed', case=False).sum()),
                     index=['induced', 'repressed'])


def get_query_result(query: str, user_lists: Optional[Tuple[pd.DataFrame, Dict]] = None,
                     cache_path: Optional[Union[str, Path]] = None) -> Tuple[pd.DataFrame, pd.DataFrame, Dict]:
    stats = {}
    result = parse_query(query)
    metadata = get_metadata(result.columns.levels[1])
    result = expand_ref_ids(result, 1)
    result = reorder_data(result)

    stats['total'] = result.loc[:, (slice(None), slice(None), slice(None), 'EDGE')].groupby(
        level=[0, 1, 2], axis=1).count().sum()

    if user_lists:
        result = result[result.index.isin(user_lists[0].index)]

    if cache_path:  # cache edges here
        result.to_pickle(cache_path + '/tabular_output.pickle.gz')
        metadata.to_pickle(cache_path + '/metadata.pickle.gz')

    result = trim_edges(result)

    counts = get_tf_count(result)

    result = pd.concat([counts, result], axis=1)

    if user_lists:
        result = user_lists[0].merge(result, left_index=True, right_index=True, how='inner')
        result = result.sort_values(['User List Count', 'User List'])
    else:
        result = pd.concat([
            pd.DataFrame(np.nan, columns=['User List', 'User List Count'], index=result.index),
            result
        ], axis=1)

    result = ANNOTATIONS.merge(result, how='right', left_index=True, right_index=True)

    result = result.sort_values('TF Count', ascending=False, kind='mergesort')

    # return statistics here as well
    return result, metadata, stats
