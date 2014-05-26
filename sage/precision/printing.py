def repr_operator(op,args):
    from sage.functions.other import floor, ceil
    l = len(args)
    if l == 0: return ""
    if l == 1: return args[0]
    terms = [ ]
    lengths = [ ]
    maxlines = 0
    for i in range(l): 
        term = args[i].split('\n') 
        length = max([ len(lgn) for lgn in term ])
        if len(term) > maxlines:   
            maxlines = len(term)
        terms.append(term)
        lengths.append(length)
    if maxlines > 1:
        if op == "": op = " "
        else:
            if op[0] != " ": op = " " + op
            if op[-1] != " ": op += " "
    shifts = [ ceil((maxlines - len(terms[i]))/2) for i in range(l) ]
    lineop = floor(maxlines/2)
    s = ""
    for n in range(maxlines):
        for i in range(l):
            m = n - shifts[i]
            if m < 0 or m >= len(terms[i]):
                t = ""
            else:
                t = terms[i][m]
            s += t + " "*(lengths[i]-len(t))
            if i < l-1:
                if n == lineop: s += op
                else: s += " "*len(op)
        s += "\n"
    return simplify(s[:-1])

def repr_exponant(s,exp):
    ls = s.split('\n')
    lexp = exp.split('\n')
    length = max([ len(line) for line in ls ])
    ns = len(ls)
    nexp = len(lexp)
    if nexp > 1: op = " "
    else: op = "^"
    s = ""
    for i in range(ns+nexp-1):
        if i < nexp-1:
            s += " " * length
        else:
            t = ls[i-nexp+1]
            s += t
            s += " " * (length - len(t))
        if i < nexp:
            s += op + lexp[i]
        s += "\n"
    return simplify(s[:-1])

def repr_index(s,ind):
    ls = s.split('\n')
    lind = ind.split('\n')
    length = max([ len(line) for line in ls ])
    ns = len(ls)
    nind = len(lind)
    if nind > 1: op = " "
    else: op = "_"
    s = ""
    for i in range(ns+nind-1):
        if i < ns:
            t = ls[i]
            s += t
            s += " " * (length - len(t))
        else:
            s += " " * length
        if i >= ns-1:
            s += op + lind[i-ns+1]
        s += "\n"
    return s[:-1]

def repr_parenthesis(s):
    lines = s.split('\n')
    l = len(lines)
    if l == 0:
        return ""
    if l == 1:
        return "(" + s + ")"
    length = max([ len(line) for line in lines ])
    sp = ""
    for n in range(l):
        lines[n] += " " * (length - len(lines[n]))
        if n == 0:
            sp += "/ " + lines[n] + " \\\n"
        elif n == l-1:
            sp += "\\ " + lines[n] + " /\n"
        else:
            sp += "| " + lines[n] + " |\n"
    return sp[:-1]

def repr_function(name,args):
    from sage.functions.other import floor
    s = repr_parenthesis(repr_operator(", ",args))
    lines = s.split('\n')
    lineop = floor(len(lines)/2)
    sp = ""
    if len(lines) > 1:
        name += " "
    for n in range(len(lines)):
        if n == lineop: 
            sp += name
        else:
            sp += " " * len(name)
        sp += lines[n] + "\n"
    return sp[:-1]

def simplify(s):
    import re
    lines = s.split("\n")
    for i in range(len(lines)):
        line = lines[i]
        positions = [ ]
        replacement = [ ]
        for m in re.finditer("(?:\+ *- ?|\* *1 */)",line):
            positions.append((m.start(),m.end()))
            if line[m.start()] == '+':
                replacement.append("- ")
            else:
                replacement.append("/")
        for n in range(len(positions)):
            (start, end) = positions[n]
            erase = True
            spaces = " " * (end-start)
            for j in range(len(lines)):
                if i == j: continue
                if lines[j][start:end] != spaces:
                    erase = False
                    break
            if not erase: break
            repl = replacement[n]
            for j in range(len(lines)):
                if i == j:
                    lines[j] = lines[j][:start] + repl + lines[j][end:]
                else:
                    lines[j] = lines[j][:start] + (" " * len(repl)) + lines[j][end:]
    s = ""
    for i in range(len(lines)):
        s += lines[i] + "\n"
    return s[:-1]
