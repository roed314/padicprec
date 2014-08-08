def my_quo_rem(self, right):
    """
    An implementation of quo_rem using lists of coefficients.

    Faster than the one provided by the class.
    """
    if right.is_zero():
        raise ZeroDivisionError, "cannot divide by a polynomial indistinguishable from 0"
    a = self.list(); da = len(a)-1
    b = right.list(); db = right.degree()
    inv = ~b[db]
    q = [ ]
    for i in range(da,db-1,-1):
        c = inv*a[i]
        q.append(c)
        for j in range(db):
            a[j+i-db] -= c*b[j]
    q.reverse()
    K = self.base_ring().fraction_field()
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
    parent = PolynomialRing(K, name=self.parent().variable_name())
    return parent(q), parent(a[:db])

def my_rem(self, right):
    _, r = my_quo_rem(self, right)
    return r

def my_quo(self, right):
    q, _ = my_quo_rem(self, right)
    return q


def val_poly(P):
    """
    Returns the valuation of the polynomial P
    """
    val = Infinity
    for c in P.list():
        if c.is_zero(): continue
        v = c.valuation()
        if v < val: val = v
    return val

def prec_poly(P):
    """
    Returns the precision of the polynomial P
    """
    prec = Infinity
    for c in P.list():
        p = c.precision_absolute()
        if p < prec: prec = p
    return prec


def xgcd_step(A,B):
    """
    Use Euclide algorithm to compute 
      - the gcd D of A, B
      - a 2x2 matrix M such that:
          [A,B] * M^t = [D,0]
    """
    S = A.parent()
    R1 = A; U1 = S(1); V1 = S(0)
    R2 = B; U2 = S(0); V2 = S(1)
    while R2 != 0:
        Q, R = my_quo_rem(R1,R2)
        U = U1 - Q*U2
        V = V1 - Q*V2
        R1, R2 = R2, R
        U1, U2 = U2, U
        V1, V2 = V2, V
    M = matrix(2,2, [U1,V1,U2,V2])
    return R1, M


def _renormalize_UV(A, B, U, V):
    """
    Compute D', U' and V' such that:

        A U' + B V' = D'

    where D' is the factor of D = AU + BV of slope 0

    Be careful: A and B have to be monic
    """
    D = A*U + B*V

    coeffs = D.list()
    for i in range(len(coeffs)-1, -1, -1):
        if coeffs[i].valuation() == 0:
            deg = i
            break
    a = D.truncate(deg+1)
    b = v = D.parent()(1)
    x = my_rem(D, a)
    while(not x.is_zero()):
        a += my_rem(v*x, a)
        b, x = my_quo_rem(D,a)
        b = my_rem(b, a)
        v = my_rem(v * (2 - b*v), a)
    Dp = a

    C = my_quo(D, Dp)
    Cinv = 1; tst = C
    prec = tst.base_ring().precision_cap()
    while val_poly(tst-1) < prec:
        Cinv = my_rem(Cinv * (2 - tst), B)
        tst = my_rem(C * Cinv, B)
    Up = my_rem(Cinv * U, B)
    Vp = my_quo(Dp - A*Up, B)
    return Dp, Up, Vp

def xgcd_stable(A,B):
    """
    Stable (and slow) algorithm for computing xgcd
    """
    S = A.parent()
    K = S.base_ring()
    prec = K.precision_cap()
    k = K.residue_field()
    Sbar = PolynomialRing(k, name='xbar')
    M = MatrixSpace(K,2)(1) 
    a = A; b = B

    while b != 0:
        # We compute the gcd mod p
        alpha = Sbar(a); beta = Sbar(b)
        delta, mu = xgcd_step(alpha,beta)

        # We lift the transformation matrix
        Mbar = matrix(2,2, [ S([ K(c).lift_to_precision(prec) for c in entry.list() ]) for entry in mu.list() ])

        # We clear final terms of the gcd
        R1, U1, V1 = _renormalize_UV(a, b, Mbar[0,0], Mbar[0,1])

        # We do an extra loop in Euclide algorithm
        U2 = Mbar[1,0]; V2 = Mbar[1,1]
        R2 = a*U2 + b*V2
        Q, R = my_quo_rem(R2,R1)
        U = my_rem(U2 - Q*U1, b)
        V = my_quo(R - a*U, b)

        # We renormalize U, V, R
        val = val_poly(R)
        if val is not Infinity:
            coeff = K(1) >> val
            R *= coeff
            U *= coeff
            V *= coeff

        # We rebuild M, a, b for the next loop
        if R.leading_coefficient().valuation() > 0 and R != 0:
            R, U, V = _renormalize_UV(a, b, U, V)
        M = matrix(2,2, [U1,V1,U,V]) * M
        a = R1; b = R

    U = M[0,0]; V = M[0,1]
    U = my_rem(U,B)
    V = my_quo(a - A*U, B)
    lc = a.leading_coefficient()
    return a/lc, U/lc, V/lc


def xgcd_euclide(A,B):
    """
    Euclide algorithm for computing xgcd
    """
    D, M = xgcd_step(A,B)
    lc = D.leading_coefficient()
    return D/lc, M[0,0]/lc, M[0,1]/lc


def xgcd_test(A,B):
    """
    Computes xgcd(A,B) with:
    - Euclide algorithm, and
    - our stable algotithm
    and outputs timings and precisions.
    """
    t = cputime()

    print "Usual Euclide algorithm"
    D, U, V = xgcd_euclide(A,B)
    print "  degree of GCD: %s" % D.degree()
    prec_euclide = min(prec_poly(U), prec_poly(V))
    print "  precision: %s" % prec_euclide
    print "  time: %s" % cputime(t)

    print "Stable algorithm"
    D, U, V = xgcd_stable(A,B)
    print "  degree of GCD: %s" % D.degree()
    prec_stable = min(prec_poly(U), prec_poly(V))
    print "  precision: %s" % prec_stable
    print "  time: %s" % cputime(t)


def xgcd_tests(p=3, prec=20, degree=20, repeat=10):
    """
    Runs a bunch of tests on random inputs
    """
    R = Zp(p, prec)
    K = R.fraction_field()
    S.<x> = PolynomialRing(K)
    for _ in range(repeat):
        A = S([ R.random_element() for _ in range(degree) ] + [ R(1) ])
        B = S([ R.random_element() for _ in range(degree) ] + [ R(1) ])
        xgcd_test(A,B)
        print "--"
