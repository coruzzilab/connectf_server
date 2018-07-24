import re
from collections import deque
from functools import partial
from pathlib import Path
from typing import Any, Dict, Optional, Sequence, Tuple, Union
from uuid import uuid4

import numpy as np
import pandas as pd
import pyparsing as pp
from django.db.models import Q

from querytgdb.models import AnalysisIddata, Annotation, DAPdata, Interactions, MetaIddata, ReferenceId, Regulation

__all__ = ['get_query_result', 'expand_ref_ids', 'ANNOTATIONS']


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
    else:
        return pd.DataFrame(False, columns=df.columns, index=df.index)


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
        elif re.match(r'^dap$', key, flags=re.I):
            return df.groupby(level=[0, 1], axis=1).apply(apply_search_column, value=value, key='DAP')
        elif re.match(r'^has_column$', key, flags=re.I):
            value = value.upper()
            return df.groupby(level=[0, 1], axis=1).apply(apply_has_column, value=value)
        else:
            return query_metadata(df, key, value)
    else:
        return query


def get_tf_data(query: str) -> TargetFrame:
    """
    Get data for single TF
    :param query:
    :return:
    """
    df = TargetFrame(
        interactions.filter(db_tf_id__db_tf_agi__exact=query).values_list(
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

        dap = TargetSeries(
            DAPdata.objects.filter(db_tfid__db_tf_agi__exact=query).values_list('ath_id__agi_id', flat=True).iterator())
        df = df.merge(reg, on=['REF', 'TARGET'], how='left')
        if not dap.empty:
            df.loc[df['TARGET'].isin(dap), 'DAP'] = 'Present'
        df = (df.pivot(index='TARGET', columns='REF')
              .swaplevel(0, 1, axis=1)
              .sort_index(axis=1, level=0, sort_remaining=False)
              .dropna(axis=1, how='all'))

    df.columns = pd.MultiIndex.from_tuples((query, *c) for c in df.columns)
    df.filter_string += query

    return df


def get_all_tf(query: str) -> TargetFrame:
    """
    Get data for all TFs at once
    :param query:
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

    dap = TargetFrame(
        DAPdata.objects.values_list('db_tfid__db_tf_agi', 'ath_id__agi_id').iterator(),
        columns=['TF', 'TARGET'])
    dap['DAP'] = 'Present'
    df = (df.merge(reg, on=['REF', 'TARGET'], how='left')
          .merge(dap, on=['TF', 'TARGET'], how='left')
          .set_index(['TF', 'REF', 'TARGET'])
          .unstack(level=[0, 1])
          .reorder_levels([1, 2, 0], axis=1)
          .sort_index(axis=1, level=[0, 1], sort_remaining=False)
          .dropna(how='all', axis=1))

    if query == 'oralltf':
        df.filter_string += 'oralltf'
    elif query == 'andalltf':
        df = df[df.loc[:, (slice(None), slice(None), 'EDGE')].notna().all(axis=1)]
        df.filter_string += 'andalltf'
    else:
        raise ValueError('invalid query')

    return df


def get_suffix(prec: TargetFrame, succ: TargetFrame) -> Tuple[str, str]:
    if prec.filter_string == succ.filter_string:
        return '_' + str(uuid4()), '_' + str(uuid4())

    return ' "' + prec.filter_string + '"', ' "' + succ.filter_string + '"'


def get_tf(query: Union[pp.ParseResults, str, TargetFrame]) -> TargetFrame:
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
                    succ = get_tf(next(it))
                    succ.include = not succ.include
                    succ.filter_string = 'not ' + succ.filter_string
                    stack.append(succ)
                elif is_modifier(curr):
                    prec = get_tf(stack.pop())
                    mod = get_mod(prec, curr)
                    prec = prec[mod].dropna(how='all')

                    # filter out empty tfs
                    prec = prec.groupby(level=[0, 1], axis=1).filter(lambda x: x.notna().any(axis=None))
                    prec.filter_string += '[' + mod_to_str(curr) + ']'

                    stack.append(prec)
                else:
                    stack.append(curr)
        except StopIteration:
            return get_tf(stack.pop())
    elif isinstance(query, (TargetFrame, TargetSeries)):
        return query
    else:
        if query in ('andalltf', 'oralltf'):
            return get_all_tf(query)

        return get_tf_data(query)


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
    df = TargetFrame(refs.values_list(
        'ref_id',
        'analysis_id__analysisiddata__analysis_type',
        'analysis_id__analysisiddata__analysis_value').iterator(),
                     columns=['REF', 'KEY', 'VALUE'])

    df = pd.concat([
        df,
        TargetFrame(
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


def parse_query(query: str) -> TargetFrame:
    parse = expr.parseString(query, parseAll=True)

    result = get_tf(parse.get('query'))

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


def get_query_result(query: str, user_lists: Optional[Tuple[TargetFrame, Dict]] = None,
                     cache_path: Optional[Union[str, Path]] = None) -> Tuple[TargetFrame, TargetFrame, Dict]:
    stats = {}
    result = parse_query(query)
    metadata = get_metadata(result.columns.get_level_values(1))
    result = expand_ref_ids(result, 1)

    stats['total'] = result.loc[:, (slice(None), slice(None), slice(None), 'EDGE')].groupby(
        level=[0, 1, 2], axis=1).count().sum()

    if cache_path:
        result.to_pickle(cache_path + '/tabular_output_unfiltered.pickle.gz')

    if user_lists:
        result = result[result.index.isin(user_lists[0].index)]

    if cache_path:  # cache edges here
        result.to_pickle(cache_path + '/tabular_output.pickle.gz')
        metadata.to_pickle(cache_path + '/metadata.pickle.gz')

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

    result = ANNOTATIONS.merge(result, how='right', left_index=True, right_index=True)

    result = result.sort_values('TF Count', ascending=False, kind='mergesort')

    # return statistics here as well
    return result, metadata, stats
