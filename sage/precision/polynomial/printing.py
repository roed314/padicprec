from precision.printing import repr_operator, repr_parenthesis, repr_exponant, simplify

def repr_polynomial(coeffs, variable_name, variable_level_parenthesis=3):
    length = len(coeffs)
    if length == 0: return '0'
    terms = [ ]
    for i in range(length):
        coeff = coeffs[i]
        s = str(coeff)
        if s == '0': continue
        if coeff.parenthesis_level() < 1:
            s = repr_parenthesis(s)
        if i == 1:
            if variable_level_parenthesis < 1:
                monome = repr_parenthesis(variable_name)
            else:
                monome = variable_name
        elif i > 1:
            if variable_level_parenthesis < 2:
                monome = repr_parenthesis(variable_name)
            else:
                monome = variable_name
            monome = repr_exponant(monome, str(i))
        if i == 0:
            terms.append(s)
        elif s == '1':
            terms.append(monome)
        else:
            terms.append(repr_operator("*", [s, monome]))
    terms.reverse()
    s = repr_operator(" + ", terms)
    if s == '': s = '0'
    return simplify(s)
