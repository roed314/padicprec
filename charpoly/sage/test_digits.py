#from digits_of_precision import *

def digits_standard_basis(M):
    R = M[0,0].base_ring()
    S = PolynomialRing(R, name='X')
    C = (M - S.gen()).adjoint()   # à remplacer par le truc de Tristan

    n = M.nrows()
    vectors = [ ]
    for i in range(n):
        for j in range(n):
            vectors.append([ C[i,j][k] for k in range(n) ])

    return digits_of_precision(vectors)


def digits_legendre_diagonal(M):
    R = M[0,0].base_ring()
    S = PolynomialRing(R, name='X')
    C = (M - S.gen()).adjoint()   # à remplacer par le truc de Tristan
    diag = M.diagonal()

    n = M.nrows()
    vectors = [ ]
    for i in range(n):
        for j in range(n):
            vectors.append([ C[i,j](x) for x in diag ])

    return digits_of_precision(vectors)

def digits_legendre(M, evaluation_points):
    R = M[0,0].base_ring()
    S = PolynomialRing(R, name='X')
    C = (M - S.gen()).adjoint()   # à remplacer par le truc de Tristan

    n = M.nrows()
    vectors = [ ]
    for i in range(n):
        for j in range(n):
            vectors.append([ C[i,j](x) for x in evaluation_points ])

    return digits_of_precision(vectors)
