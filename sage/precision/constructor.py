from precision.parent_inexact import Parent_inexact


def PolynomialRing(parent,*args,**kwargs):
    if isinstance(parent,Parent_inexact):
        from precision.polynomial.polynomial_ring import PolynomialRing_inexact as constructor
    else:
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing as constructor
    return constructor(parent,*args,**kwargs)


def MatrixSpace(parent,*args,**kwargs):
    if isinstance(parent,Parent_inexact):
        from precision.matrix.matrix_space import MatrixSpace_inexact as constructor
    else:
        from sage.matrix.matrix_space import MatrixSpace as constructor
    return constructor(parent,*args,**kwargs)
