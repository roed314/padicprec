R = Zp(2,prec=1000)
prec = 500
S.<x> = PolynomialRing(R)

# Precision functions
#####################

def precmoy_absolute(P):
    prec = 0; count = 0
    for c in P.list():
        pr = c.precision_absolute()
        if pr < Infinity:
            prec += pr
            count += 1
    return prec/count

def valmoy_absolute(P):
    val = 0; count = 0
    for c in P.list():
        v = c.valuation()
        if v < Infinity:
            val += v
            count += 1
    return val/count

def precmoy_relative(P):
    prec = 0; count = 0
    for c in P.list():
        pr = c.precision_relative()
        if pr < Infinity:
            prec += pr
            count += 1
    return prec/count

# Euclidean division
####################

def my_quo_rem(self, right):
    """
    An implementation of quo_rem using lists of coefficients.

    Faster than the one provided by the class.
    """
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


# Main
######

def test(modulus, N=100, repeat=1):
    modulus = S(modulus)
    d = modulus.degree()
    sval = sabs = stheory = 0

    x = S.gen()
    Basis = [ ]
    for i in range(d):
        Basis.append(x**i)
    for _ in range(repeat):
        Lattice = [ ]
        for i in range(d):
            Lattice.append(x**i)
        res = S(1)
        theory = 0
        for c in range(N):
            A = S([ R.random_element().add_bigoh(prec) for i in range(d) ])
            Lat = [ my_rem(res*dA,modulus) for dA in Basis ] + [ my_rem(dres*A,modulus) for dres in Lattice ]
            res = my_rem(res*A, modulus)
            Lattice = [ ]
            for i in range(d):
                val = Infinity; index = -1
                for ind in range(len(Lat)):
                    if Lat[ind][i] == 0: continue
                    v = Lat[ind][i].valuation()
                    if v < val: 
                        val = v; index = ind
                if index == -1: raise RuntimeError
                gen = Lat[index]
                Lattice.append(gen)
                del Lat[index]
                if c == N-1:
                    theory += gen[i].valuation()
                for ind in range(len(Lat)):
                    scalar = Lat[ind][i] // gen[i]
                    scalar = scalar.lift_to_precision(prec)
                    Lat[ind] -= scalar * gen
        val = RR(valmoy_absolute(res))
        sval += val
        print ". Valuation: %s" % val
        abs = RR(precmoy_absolute(res)) - prec
        sabs += abs
        print ". Absolute: %s" % abs
        theory = RR(theory/d)
        stheory += theory
        print ". Theory: %s" % theory

    print "* Valuation: %s" % (sval/repeat)
    print "* Absolute: %s" % (sabs/repeat)
    print "* Theory: %s" % (stheory/repeat)
