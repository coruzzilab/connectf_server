import pyparsing as pp

__all__ = ['meta_parser']

eq = pp.Literal('=')
ident = pp.Word(pp.alphanums + '_') + pp.FollowedBy(eq)
p_ident = pp.CaselessLiteral('P-value') + pp.FollowedBy(eq)
value = pp.Word(pp.alphanums + '_-.')
oper = pp.locatedExpr(pp.oneOf('And Or AndNot', caseless=True))("oper*")

criteria = pp.Group(ident + pp.Suppress(eq) + value)
p_criteria = pp.locatedExpr(pp.Group(p_ident + pp.Suppress(eq) + value))('p_value*')

meta_parser: pp.ParserElement = pp.infixNotation(criteria ^ p_criteria,
                                                 [(oper, 2, pp.opAssoc.LEFT)],
                                                 lpar=pp.Suppress('['),
                                                 rpar=pp.Suppress(']'))
