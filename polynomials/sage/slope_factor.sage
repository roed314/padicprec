load("mulmod.sage")
p = 11
R = Zp(11,prec=100)
prec = 100
SZ.<X> = PolynomialRing(ZZ)
S.<x> = PolynomialRing(R)

########################################
#     Tools
########################################

def prec_absolute(P):
    prec = 0; count = 0
    for c in P.list():
        pr = c.precision_absolute()
        if pr < Infinity:
            prec += pr
            count += 1
    return prec

########################################
#Here for the slope factorisation
########################################

def slope_factor(P,d,niter,S):
    """
	Compute the first factor of the factorization of P into two polynomials along 
        the vertex at d of the Newton Polygon of P
    """
    x=S.gen()
    A0 = sum([P.list()[i]*x^i for i in range(d+1)])
    A0=A0*A0.list()[d]^(-1)
    Ai=A0
    V=S(1)
    for i in range(niter):
        (Q,R)=my_quo_rem(V*P,Ai)
        Ai=(Ai+R).truncate(d+1)
        (Q,R)=my_quo_rem(P,Ai)
        V=my_rem(2*V-V^2*Q,Ai)
    return Ai,V

########################################
# Design of test of slope factorization
########################################

def random_poly_above_NP(NP,R, truncprec):
    """
	create a polynomial with random coefficients 
        above the Newton polygon     NP.
        highest coefficient is of the form p^l
    """
    p=R.uniformizer()
    vert = NP.vertices()
    d = vert[len(vert)-1][0]
    P = p^(ceil(NP(0)))
    for i in range(1,d):
        P += ((R.random_element()*p^(ceil(NP(i)))).add_bigoh(truncprec))*x^i
    P+=(R(p^(ceil(NP(d)))).add_bigoh(truncprec))*x^d
    return P
	
def test_factor(NP1,NP2,S,niter,truncprec):
    """
        create two random polynomials A and B with Newton polygons above NP1 and NP2
        Compute the product AB and factor it by slope factorisation to obtain A
        Compare precision on A, A2 and 
    """
    R=S.base_ring()
    x=S.gen()
    A = random_poly_above_NP(NP1,R,truncprec)
    B = random_poly_above_NP(NP2,R,truncprec)
    P = A*B

    vert = NP1.vertices()
    d = vert[len(vert)-1][0]
    degB = NP2.vertices()[len(NP2.vertices())-1][0]
    A2,V = slope_factor(P,d,niter,S)

    A0 = sum([P.list()[i]*x^i for i in range(d+1)])

    ##################
    # Jagged precision
    ##################
    abs = prec_absolute(A2) - (d+1)*truncprec
    print("jagged precision")
    print abs

    ##################
    # Newton precision
    ##################
    precA = NewtonPolygon([(0,truncprec),(d-1,truncprec)])
    precB = NewtonPolygon([(0,truncprec),(degB-1,truncprec)])
    precP = precNewton_mul(A,precA,B,precB)
    #precres = precNewton_mul(P,precP,V,precision_polygon(V))
    #precres = precNewton_mul(P,precP,V,precision_polygon(V))
    #  Using Psi   
    precres = precP*V.newton_polygon()
    quo, precres = Newton_quorem(precres,A0.newton_polygon())
    nwton = sum([ceil(precres(i)) for i in range(d)]) - d*truncprec
    print("Newton precision")
    print nwton

    ####################    
    # Lattice estimation
    ####################    
    dA = [x^i for i in range(d+1)]
    dB = [x^i for i in range(NP1.vertices()[len(NP1.vertices())-1][0]+1)]
    dP = [ B*da for da in dA] + [ db*A for db in dB ]
    dA2 = [ my_rem(V*dp,A) for dp in dP]

    # SNF of dA2
    theory = theory2 = 0
    Lat = dA2
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
        theory += gen[i].valuation()
        theory2 += min([elt[i].valuation() for elt in Lattice])
        for ind in range(len(Lat)):
            scalar = Lat[ind][i] // gen[i]
            scalar = scalar.lift_to_precision(truncprec)
            Lat[ind] -= scalar * gen
    print("lattice precision")
    print theory, theory2, theory-theory2


    return A,B,A2,V


########################################
# Actual testing
########################################


#NPtest1 =NewtonPolygon([(0,2),(4,0)])
#NPtest2 =NewtonPolygon([(0,0),(4,3)])

NPtest1 =NewtonPolygon([(0,5),(2,-2),(6,0)])
NPtest2 =NewtonPolygon([(0,0),(4,3),(7,10)])

A,B,A2,V=test_factor(NPtest1,NPtest2,S,50, 100)
P=A*B
#nana=(A-A2).newton_polygon().plot()
Q,R=my_quo_rem(P,A)
Q2,R2=my_quo_rem(P,A2)

######

NPtest1 =NewtonPolygon([(0,5),(2,1),(6,0)])
NPtest2 =NewtonPolygon([(0,0),(3,0),(5,10)])

A,B,A2,V=test_factor(NPtest1,NPtest2,S,50, 100)


#####
NPtest1 =NewtonPolygon([(0,5),(2,1),(3,0)])
NPtest2 =NewtonPolygon([(0,0),(5,-2),(7,10)])

A,B,A2,V=test_factor(NPtest1,NPtest2,S,50, 100)


