def val_poly(P):
    val = Infinity
    for c in P.list():
        if c.is_zero(): continue
        v = c.valuation()
        if v < val: val = v
    return val

def prec_poly(P):
    prec = Infinity
    for c in P.list():
        p = c.precision_absolute()
        if p < prec: prec = p
    return prec

def lift_poly(P, prec):
    parent = P.parent()
    return parent([ c.lift_to_precision(prec) for c in P.list() ])

def add_bigoh_poly(P, prec):
    parent = P.parent()
    return parent([ c.add_bigoh(prec) for c in P.list() ])

def my_rem_integral(self, right):
    """
    Compute the remainder in the Euclidean division
    of:

        lc(right)^(deg(self) - deg(right)) * self

    by right.

    (It avoids all divisions.)
    """
    a = self.list(); da = len(a)-1
    b = right.list(); db = right.degree()
    lb = b[db]
    for i in range(da,db-1,-1):
        c = a[i]
        for j in range(i-db):
            a[j] *= lb
        for j in range(db):
            a[j+i-db] = a[j+i-db]*lb - c*b[j]
    return self.parent()(a[:db])


def subresultants_stable(A, B):
    """
    Returns all subresultants of A and B (until one of them vanishes).

    NB: For now, A and B must be monic, have integral coefficients 
    and share the same degree.
    """
    d = A.degree()
    if B.degree() != d:
        raise NotImplementedError("A and B does not have same degree")
    if A.leading_coefficient() != 1 or B.leading_coefficient() != 1:
        raise NotImplementedError("A or B is not monic")
    if val_poly(A) < 0 or val_poly(B) < 0:
        raise NotImplementedError("A or B does not have integral coefficients")
    prec = min(prec_poly(A), prec_poly(B))
    R0 = A; R1 = my_rem_integral(A,B)
    sres = [ ]
    while val_poly(R1) < prec:
        sres.append(R1)
        d -= 1
        if R1[d] == 0:
            raise NotImplementedError
        a0 = R0.leading_coefficient(); v0 = a0.valuation()
        a1 = R1.leading_coefficient(); v1 = a1.valuation()
        R0 = lift_poly(R0, prec + 2*v0 + 2*v1)
        R1 = lift_poly(R1, prec + 2*v0 + 2*v1)
        R = my_rem_integral(R0, R1)
        R *= ~(R0.leading_coefficient())^2
        R0, R1 = R1, R
    return [ add_bigoh_poly(R, prec) for R in sres ]


def gcd_stable(A,B):
    """
    Stable algorithm to compute gcd(A,B)
    """
    D = subresultants_stable(A,B)[-1]
    return D.monic()


# Tests

def test(p=3, prec=20, degree=20):
    R = Zp(p, prec + 10*ceil(log(degree)/log(p)))

    S.<x> = R[]
    A = S([ R.random_element().add_bigoh(prec) for _ in range(degree) ] + [ R(1).add_bigoh(prec) ])
    B = S([ R.random_element().add_bigoh(prec) for _ in range(degree) ] + [ R(1).add_bigoh(prec) ])

    print "A ="; print A; print
    print "B ="; print B; print

    print "gcd_stable(A*A, A*B) ="
    t = cputime()
    D = gcd_stable(A*A, A*B)
    t = cputime(t)
    print D
    print "Time: %s" % t
    print

    print "gcd(A*A, A*B) ="
    t = cputime()
    D, _, _ = xgcd(A*A, A*B)
    t = cputime(t)
    print D
    print "(Degree = %s)" % D.degree()
    print "Time: %s" % t


# Some statistics

def stat_subresultants(p=3, prec=20, degree=20, repeat=10):
    R = Zp(p, prec + 10*ceil(log(degree) / log(p)))
    S.<x> = PolynomialRing(R)
    counts = [ ]; count = 0
    for _ in range(repeat):
        A = S([ R.random_element().add_bigoh(prec) for _ in range(degree) ] + [ R(1).add_bigoh(prec) ])
        B = S([ R.random_element().add_bigoh(prec) for _ in range(degree) ] + [ R(1).add_bigoh(prec) ])
        for P in subresultants_stable(A,B):
            d = P.degree()
            for _ in range(len(counts), d+1):
                counts.append([])
            v = P.leading_coefficient().valuation()
            for _ in range(len(counts[d]), v+1):
                counts[d].append(0)
            counts[d][v] += 1
        count += 1
    print "EXPECTED VALUES FOR THE VALUATION"
    sesp = 0
    for d in range(len(counts)):
        esp = 0
        for v in range(len(counts[d])):
            esp += v*counts[d][v]
        print "Degree %s: %s" % (d, RR(esp/count))
        sesp += esp
    print "Sum: %s" % RR(sesp/count)



# Under construction

def _subresultants(A, B, prec_lift=None):
    prec = min(prec_poly(A), prec_poly(B))
    if prec_lift is None:
        prec_lift = A.base_ring().precision_cap()
    if A.degree() < B.degree():
        R = [ lift_poly(B, prec_lift), lift_poly(A, prec_lift) ]
    else:
        R = [ lift_poly(A, prec_lift), lift_poly(B, prec_lift) ]
    d = [ "", R[0].degree() - R[1].degree() ]
    a = [ 1, R[1].leading_coefficient() ]
    c = [ 1, a[1]^d[1] ]
    i = 1
    minprec = Infinity
    MS = MatrixSpace(A.parent().change_ring(A.base_ring().fraction_field()),2,2)
    M = MS(1)
    while R[i] != 0:
        nQ, nR = my_quo_rem(a[i]^(d[i]+1)*R[i-1], R[i]) 
        coeff = ~(a[i-1] * c[i-1]^d[i])
        nQ *= coeff; nR *= coeff
        minprec = min(minprec, prec_poly(nR))
        if minprec < 0:
            raise RuntimeError
        R.append(lift_poly(nR, prec_lift))
        M = MS([0, 1, (a[i]*c[i])/(a[i-1]*c[i-1]), nQ]) * M
        print [ val_poly(entry) for entry in M.list() ], val_poly(M.determinant())
        i += 1
        d.append(R[i-1].degree() - R[i].degree())
        a.append(R[i].leading_coefficient())
        c.append(a[i]^d[i] / c[i-1]^(d[i]-1))
    print R[-2]
    return minprec, [ add_bigoh_poly(P,prec) for P in R[2:-1] ]
