def repr_matrix(coeffs):
    from sage.functions.other import floor
    ncols = 0
    lengths = [ ]
    for row in coeffs:
        l = len(row)
        if ncols < l:
            lengths.extend([0] * (l-ncols))
            ncols = l
        for j in range(l):
            elt = str(row[j]); l = len(elt)
            if lengths[j] < len(elt)+2:
                lengths[j] = len(elt)+2
    s = ""
    for row in coeffs:
        s += "["
        for j in range(len(row)):
            elt = str(row[j]); l = len(elt)
            nbspace = floor((lengths[j]-l)/2)
            s += " " * nbspace
            s += elt
            s += " " * (lengths[j]-l-nbspace)
        for j in range(len(row),ncols):
            s += " " * lengths[j]
        s += "]\n"
    return s[:-1]

