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
    R0 = A; R1 = -my_rem_integral(A,B)
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


def discriminant_stable(P):
    Q = P.derivative() + P
    disc = subresultants_stable(P,Q)[-1]
    if disc.degree() > 0:
        return P.parent().base_ring()(0)
    else:
        return disc[0]


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

def twisted_convolution(jlaw, q):
    L = len(jlaw)
    res = [ ]
    for s in range(L):
        res_s = [ ]
        for x in range(s+1):
            prob = 0
            if x <= s/2:
                prob = q^(-x) * jlaw[s-x][x]
                for xp in range(x+1, s-x+1):
                    prob += (1 - 1/q) * q^(-x) * jlaw[s-x][xp]
            res_s.append(prob)
        res.append(res_s)
    return res

def compute_law(L, p=2, max=10):
    L = list(L)
    L.sort()
    jlaw = None
    for j in range(len(L)):
        if jlaw is None:
            di = L[j]
            jlaw = [ n*[0] + [ (1 - 1/p^di) * p^(-di*n) ] for n in range(max) ]
        else:
            di = L[j] - L[j-1]
            jlaw = twisted_convolution(jlaw, p^di)
    return jlaw # [ sum(l) for l in jlaw ]


def stat_subresultants(p=3, prec=20, degree=20, repeat=1000):
    R = Zp(p, 5*prec)
    S.<x> = PolynomialRing(R)
    counts = [ ]; count = 0
    for _ in range(degree): counts.append([ ])
    for _ in range(repeat):
        A = S([ R.random_element().add_bigoh(prec) for _ in range(degree) ] + [ R(1).add_bigoh(prec) ])
        B = S([ R.random_element().add_bigoh(prec) for _ in range(degree) ] + [ R(1).add_bigoh(prec) ])
        j = 0
        for P in subresultants_stable(A,B):
            v = P.leading_coefficient().valuation()
            for _ in range(len(counts[j]), v+1):
                counts[j].append(0)
            counts[j][v] += 1
            j += 1
        count += 1
        if (count % 10000 == 0): print count

    odd = range(1, degree+1, 2)
    even = range(2, degree+1, 2); even.reverse()
    perm = odd + even
    expec_conj = RR(0)

    quantiles = [ None, 3.481, 5.991, 7.815, 9.488, 11.07, 12.592, 14.067, 15.507, 16.919, 18.307, 
                  19.675, 21.026, 22.362, 23.685, 24.996, 26.296, 27.587, 28.869, 30.144, 31.41, 
                  32.671, 33.924, 35.172, 36.415, 37.652, 38.885, 40.113, 41.337, 42.557, 43.773 ]

    for j in range(degree):
        print "Variable X_%s (degree %s)" % (j, degree-1-j)
        stats = [ ]; s = 0
        for c in counts[j]:
            if c < 10: break
            stats.append(c/count)
            s += c
        c = count - s
        if c < 10:
            stats[-1] += c/count
        else:
            stats.append(c/count)
        l = len(stats)
        law = compute_law(perm[:j+1], p=p, max=l-1)
        law.append(1 - sum(law))
        L = 0
        for i in range(l):
            if i < l-1:
                print " . Value %s: %s (conj: %s)" % (i, RR(stats[i]), RR(law[i]))
            else:
                print " . Value >%s: %s (conj: %s)" % (i-1, RR(stats[i]), RR(law[i]))
            L += count * (stats[i] - law[i])^2 / law[i]
        print " . L = %s" % RR(L)
        print " . quantile: %s" % quantiles[l-1]

    return counts


    #print "EXPECTED VALUES FOR THE VALUATION"
    #for j in range(degree):
    #    expec = 0; var = 0
    #    for v in range(len(counts[j])):
    #        expec += v * counts[j][v]
    #        var += v^2 * counts[j][v]
    #    expec /= count; var /= count
    #    var -= expec^2
    #    expec_conj += 1/(p^perm[j] - 1)
    #    deviation = RR(sqrt(var))
    #    print "Variable X_%s (degree %s)" % (j, degree-1-j)
    #    print " . expected value: %s (conj: %s)" % (RR(expec), expec_conj)
    #    print " . standard deviation: %s" % deviation
    #    T = RR(sqrt(repeat)) * (expec - expec_conj) / deviation
    #    print " . |T| = %s" % T.abs()


def stat_jointlaw(p=3, prec=20, degree=20, repeat=1000, cutoff=20, js=None):
    R = Zp(p, 5*prec)
    S.<x> = PolynomialRing(R)
    counts = [ ]; count = 0
    if js is None: js = range(degree)
    counts = { }
    for count in range(repeat):
        A = S([ R.random_element().add_bigoh(prec) for _ in range(degree) ] + [ R(1).add_bigoh(prec) ])
        B = S([ R.random_element().add_bigoh(prec) for _ in range(degree) ] + [ R(1).add_bigoh(prec) ])
        res = subresultants_stable(A,B)
        key = tuple([ res[degree-1-j].leading_coefficient().valuation() for j in js ])
        if counts.has_key(key):
           counts[key] += 1
        else:
           counts[key] = 1
        if (count % 10000 == 0): print count

    for t in counts.keys():
        if counts[t] < cutoff: del counts[t]
    return counts

def stat_valuation(p=3, prec=20, degree=20, repeat=1):
    R = Zp(p, 5*prec)
    S.<x> = PolynomialRing(R)
    for count in range(repeat):
        A = S([ R.random_element().add_bigoh(prec) for _ in range(degree) ] + [ R(1).add_bigoh(prec) ])
        B = S([ R.random_element().add_bigoh(prec) for _ in range(degree) ] + [ R(1).add_bigoh(prec) ])
        res = subresultants_stable(A,B)
        for P in res:
            print P.leading_coefficient().valuation(), val_poly(P)

def stat_discriminant(p=2, prec=20, degree=20, repeat=100):
    R = Zp(p, 5*prec)
    S.<x> = PolynomialRing(R)
    law = { }
    for count in range(repeat):
        P = S([ R.random_element().add_bigoh(prec) for _ in range(degree) ] + [ R(1).add_bigoh(prec) ])
        disc = discriminant_stable(P)
        val = disc.valuation()
        if law.has_key(val):
            law[val] += 1
        else:
            law[val] = 1
    return law


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


def my_quorem_integral(self, right):
    a = self.list(); da = len(a)-1
    b = right.list(); db = right.degree()
    lb = b[db]
    v = [ ]
    for i in range(da,db-1,-1):
        c = a[i]
        for j in range(i-db):
            a[j] *= lb
        for j in range(db):
            a[j+i-db] = a[j+i-db]*lb - c*b[j]
        v = [ lb*x for x in v ] + [ -c ]
    v.reverse()
    return self.parent()(a[:db]), lb^(da-db+1), self.parent()(v)

def subresultants(A, B):
    d = A.degree()
    R0 = A; 
    if B.degree() == d and (B-A).degree() < d:
        R1 = A.parent()((B-A).list()[:d])
        U0 = 1; U1 = -1
    elif B.degree() > d:
        raise NotImplementedError
    else:
        R1 = B
        U0 = 1; U1 = 0
    sres = [ ]
    cof = [ ]
    first = True
    while R1 != 0:
        sres.append(R1)
        cof.append(U1)
        d -= 1
        if R1[d] == 0:
            raise NotImplementedError
        R, U, V = my_quorem_integral(R0, R1)
        if first:
            scalar = 1
            first = False
        else:
            scalar = ~(R0.leading_coefficient())^2
        U0, U1 = U1, scalar*(U*U0 + V*U1)
        R0, R1 = R1, scalar*R
    sres.reverse(); cof.reverse()
    return sres, cof
