def val_vector(x):
    return min([ c.valuation() for c in x.list() ])
def prec_vector(x):
    return min([ c.precision_absolute() for c in x.list() ])

def stat_vectorspace(p=2, prec=100, iter=50, repeat=1000):
    K = Qp(p, prec=3*prec)
    Kint = K.integer_ring()
    V = K^3; Vint = Kint^3

    HOM = End(V)
    MS = MatrixSpace(Kint,3)

    theory = [ ]; theory2 = [ ]; practice = [ ]
    for no in range(repeat):
        L = VectorSpace(K,3)
        zero = K(0, 2*prec)
        one = K(1, 2*prec)
        delta = K(2^prec, 2*prec)
        L._add_in_place(V([one,zero,zero]))
        L1 = VectorSpace(K,3)
        L1._add_in_place(V([one,delta,zero]))
        L2 = VectorSpace(K,3)
        L2._add_in_place(V([one,zero,delta]))
        for i in range(iter):

            f1 = HOM(MS.random_element())
            f2 = HOM(MS.random_element())
            f3 = HOM(MS.random_element())
            f4 = HOM(MS.random_element())

            H1 = L.apply(f1) + L.apply(f2)
            H2 = L.apply(f3) + L.apply(f4)
            L = H1.intersection(H2)
            if L.dimension() != 1: break

            H1 = L1.apply(f1) + L1.apply(f2)
            H2 = L1.apply(f3) + L1.apply(f4)
            L1 = H1.intersection(H2)

            H1 = L2.apply(f1) + L2.apply(f2)
            H2 = L2.apply(f3) + L2.apply(f4)
            L2 = H1.intersection(H2)

        indices = L._indices
        v = L._gens[0]
        v1 = v - L1._gens[0]
        v2 = v - L2._gens[0]

        M = matrix(2,2,[v1[indices[1]], v1[indices[2]], v2[indices[1]], v2[indices[2]]])
        print [ x.valuation() for x in M.list() ]
        print RR(M.determinant().valuation()/2)
        theory.append( floor(M.determinant().valuation()/2) - prec )
        theory2.append( min(v1[indices[1]].valuation(), v2[indices[1]].valuation()) - prec )
        practice.append( prec_vector(v) - 2*prec )

        if no % 50 == 0: print no

    pts_theory = [ ]
    for i in range(min(theory), max(theory)+1):
        pts_theory.append((i, theory.count(i)))
    pts_theory2 = [ ]
    for i in range(min(theory2), max(theory2)+1):
        pts_theory2.append((i, theory2.count(i)))
    pts_practice = [ ]
    for i in range(min(practice), max(practice)+1):
        pts_practice.append((i, practice.count(i)))
    return pts_theory, pts_theory2, pts_practice
