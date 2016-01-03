import io, sys
from sage.geometry.newton_polygon import NewtonPolygon

R = Zp(2,prec=1000)
prec = 500
SZ.<X> = PolynomialRing(ZZ)
S.<x> = PolynomialRing(R)

# Newton polygons
#################

def Legendre_transform(NP, slope):
    return max([ x*slope-y for (x,y) in NP.vertices() ])

def decpolygon(self,shift):
    vertices = [ (x-shift,y) for (x,y) in self.vertices() ]
    return NewtonPolygon(vertices)

def subpolygon(self,first=0,last=Infinity):
    vertices = [ (x,y) for (x,y) in self.vertices() if x > first and x < last ] + [ (first,self(first)) ]
    if last is not Infinity:
        vertices += [ (last,self(last)) ]
    return NewtonPolygon(vertices)

def Newton_quorem(self,other):
    slopes = other.slopes(repetition=False)
    slope = slopes[-1]
    vertex = other.vertices()[-1]
    degree = vertex[0]
    Amin = subpolygon(self, last=degree-1)
    Amax = subpolygon(self, first=degree)
    if len(Amax.vertices()) == 0:
        return NewtonPolygon(), self
    shift = Legendre_transform(Amax, slope) + vertex[1] - slope*degree
    polygon = Amax + (other >> shift)
    quo = decpolygon(subpolygon(polygon, first=degree), degree)
    rem = Amin + subpolygon(polygon, last=degree-1)
    return quo,rem


# Precision functions
#####################

def precision_polygon(P):
    return NewtonPolygon([ (i,P[i].precision_absolute()) for i in srange(P.degree()+1) ])

def prec_absolute(P):
    prec = 0; count = 0
    for c in P.list():
        pr = c.precision_absolute()
        if pr < Infinity:
            prec += pr
            count += 1
    return prec

def val_absolute(P):
    val = 0; count = 0
    for c in P.list():
        v = c.valuation()
        if v < Infinity:
            val += v
            count += 1
    return val

def prec_relative(P):
    prec = 0; count = 0
    for c in P.list():
        pr = c.precision_relative()
        if pr < Infinity:
            prec += pr
            count += 1
    return prec

def precNewton_mul(f,precf,g,precg):
    return f.newton_polygon()*precg + precf*g.newton_polygon()

def precNewton_quorem(a,preca,b,precb,quo=None):
    NPa = a.newton_polygon()
    NPb = b.newton_polygon()
    if not(preca >= NPa and precb >= NPb):
        raise RuntimeError("Newton precision: not enough precision")
    if quo is None:
        quo = my_quo(a,b)
    prec = preca + precb*quo.newton_polygon()
    return Newton_quorem(prec,NPb)

def precNewton_rem(a,preca,modulus):   # the precision on modulus is infinite
    NPa = a.newton_polygon()
    if not(preca >= NPa):
        raise RuntimeError("Newton precision: not enough precision")
    quo,rem = Newton_quorem(preca,modulus.newton_polygon())
    return rem

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

def test(Modulus, N=100, repeat=1, fh=sys.stdout):
    modulus = S(Modulus)
    d = modulus.degree()
    sval = sabs = snwton = stheory = stheory2 = 0

    x = S.gen()
    Basis = [ ]
    for i in range(d):
        Basis.append(x**i)
    for no in range(repeat):
        print "No %s:" % (no+1)
        Lattice = [ ]
        for i in range(d):
            Lattice.append(x**i)
        res = S(1)
        theory = theory2 = 0
        precNewton = NewtonPolygon()
        precA = NewtonPolygon([(0,prec),(d-1,prec)])
        for c in range(N):
            A = S([ R.random_element().add_bigoh(prec) for i in range(d) ])
            precNewton = precNewton_mul(A,precA,res,precNewton)
            Lat = [ my_rem(res*dA,modulus) for dA in Basis ] + [ my_rem(dres*A,modulus) for dres in Lattice ]
            res = res*A
            precNewton = precNewton_rem(res,precNewton,modulus)
            res = my_rem(res,modulus)
            Lattice = [ ]
            for i in range(d):
                val = Infinity; index = -1
                for ind in range(len(Lat)):
                    if Lat[ind][i] == 0: continue
                    v = Lat[ind][i].valuation()
                    if v < val: 
                        val = v; index = ind
                if index == -1: raise RuntimeError("Lattice precision")
                gen = Lat[index]
                Lattice.append(gen)
                del Lat[index]
                if c == N-1:
                    theory += gen[i].valuation()
                    theory2 += min([elt[i].valuation() for elt in Lattice])
                for ind in range(len(Lat)):
                    scalar = Lat[ind][i] // gen[i]
                    scalar = scalar.lift_to_precision(prec)
                    Lat[ind] -= scalar * gen
        val = val_absolute(res)
        sval += val
        print ". Valuation: %s" % val
        abs = prec_absolute(res) - d*prec
        sabs += abs
        print ". Jagged: %s" % [ c.precision_absolute() - prec for c in res.list() ]
        nwton = sum([ceil(precNewton(i)) for i in range(d)]) - d*prec
        snwton += nwton
        print ". Newton: %s" % [ ceil(precNewton(i))-prec for i in range(d) ]
        stheory += theory
        stheory2 += theory2
        print ". Lattice: %s (%s + %s)" % (theory, theory2, theory-theory2)

    average = "* n = %s\n" % N
    average += "* modulus = %s\n" % Modulus
    average += "* Size: %s\n" % repeat
    average += "* Valuation: %s\n" % RR(sval/repeat)
    average += "* Jagged: %s\n" % RR(sabs/repeat)
    average += "* Newton: %s\n" % RR(snwton/repeat)
    average += "* Lattice: %s (%s + %s)\n\n" % (RR(stheory/repeat), RR(stheory2/repeat), RR((stheory-stheory2)/repeat))
    print "Average:"
    print average
    fh.write(unicode(average))
    fh.flush()


fh = io.open("results.txt", "a")
for modulus in [ X^5+1, X^5+X^2+1, X^5+2, (X+1)^5+2, X^5+X+2 ]:
    for N in [ 10, 50, 100 ]:
        test(modulus, N, 500, fh)
