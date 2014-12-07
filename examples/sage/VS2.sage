def val_vector(x):
    return min([ c.valuation() for c in x.list() ])
def prec_vector(x):
    return min([ c.precision_absolute() for c in x.list() ])

def stat_vectorspace(p=2, prec=100, iter=50, repeat=1000):
    K = Qp(p, prec=2*prec)
    Kint = K.integer_ring()
    V = K^3; Vint = Kint^3

    HOM = End(V)
    MS = MatrixSpace(Kint,3)

    theory = practice = 0
    for no in range(repeat):
        f1 = HOM(MS.random_element())
        f2 = HOM(MS.random_element())
        f3 = HOM(MS.random_element())
        f4 = HOM(MS.random_element())

        L = VectorSpace(K,3)
        L._add_in_place(V([1,0,0]))
        for i in range(iter):
            H1 = L.apply(f1) + L.apply(f2)
            H2 = L.apply(f3) + L.apply(f4)
            L = H1.intersection(H2)
            if L.dimension() != 1: break
        if L.dimension() != 1:
            print "Iteration %s: Not enough precision" % no
            break
        v = L._gens[0]

        L = VectorSpace(K,3)
        L._add_in_place(V([1,2^prec,0]))
        for i in range(iter):
            H1 = L.apply(f1) + L.apply(f2)
            H2 = L.apply(f3) + L.apply(f4)
            L = H1.intersection(H2)
            if L.dimension() != 1: break
        v1 = v - L._gens[0]

        L = VectorSpace(K,3)
        L._add_in_place(V([1,0,2^prec]))
        for i in range(iter):
            H1 = L.apply(f1) + L.apply(f2)
            H2 = L.apply(f3) + L.apply(f4)
            L = H1.intersection(H2)
            if L.dimension() != 1: break
        v2 = v - L._gens[0]

        theory += min(val_vector(v1), val_vector(v2)) - prec
        practice += prec_vector(v) - 2*prec

    return RR(theory/repeat), RR(practice/repeat)
