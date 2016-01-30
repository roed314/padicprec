load("mulmod.sage")
p = 2
R = Zp(p,prec=100)
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

#truncprec is not None means we are forcing flat precision on A and B

def random_poly_above_NP(NP,R,S, truncprec = None):
    """
	create a polynomial with random coefficients 
        above the Newton polygon     NP.
        highest coefficient is of the form p^l
    """
    x = S.gen()
    p=R.uniformizer()
    vert = NP.vertices()
    d = vert[len(vert)-1][0]
    prectab = []
    P = S(0)
    for i in range(1,d):
        if truncprec is None:
            c =  ((R.random_element()*p^(ceil(NP(i)))))
            prectab.append(c.precision_absolute())
            P +=c*x^i
        else:
            c =  ((R.random_element()*p^(ceil(NP(i))))).add_bigoh(truncprec)
            prectab.append(truncprec)
            P +=c*x^i              
    if truncprec is not(None):
        P+=(R(p^(ceil(NP(d)))).add_bigoh(truncprec))*x^d+ R(p^(ceil(NP(0)))).add_bigoh(truncprec)
        prectab.append(truncprec)
    else:
        P+=(R(p^(ceil(NP(d)))))*x^d+p^(ceil(NP(0)))
        prectab.append(P.list()[0].precision_absolute())
    return P,sum(prectab)
	
def test_factor(NP1,NP2,S,niter,truncprec = None):
    """
        create two random polynomials A and B with Newton polygons above NP1 and NP2
        Compute the product AB and factor it by slope factorisation to obtain A
        Compare precision on A, A2 and 
    """
    R=S.base_ring()
    x=S.gen()
    A,absA = random_poly_above_NP(NP1,R,S,truncprec)
    B, absB = random_poly_above_NP(NP2,R,S,truncprec)
    P = A*B

    vert = NP1.vertices()
    d = vert[len(vert)-1][0]
    degB = NP2.vertices()[len(NP2.vertices())-1][0]
    A2,V = slope_factor(P,d,niter,S)


    A0 = sum([P.list()[i]*x^i for i in range(d+1)])



    ##################
    # Jagged precision
    ##################
    abs = prec_absolute(A2)-A2.list()[d].precision_absolute() - absA



    ##################
    # Newton precision
    ##################
    if truncprec is not(None):
        precA = NewtonPolygon([(0,truncprec),(d-1,truncprec)])
        precB = NewtonPolygon([(0,truncprec),(degB-1,truncprec)])
        precP = precNewton_mul(A,precA,B,precB)
        #precres = precNewton_mul(P,precP,V,precision_polygon(V))
        #precres = precNewton_mul(P,precP,V,precision_polygon(V))
        #  Using Psi   
        precres = precP*V.newton_polygon()
        quo, precres = Newton_quorem(precres,A0.newton_polygon())
        nwton = sum([ceil(precres(i)) for i in range(d)]) - absA
    else:
        precA = NewtonPolygon([(i,A.list()[i].precision_absolute()) for i in range(d)])
        precB = NewtonPolygon([(i,B.list()[i].precision_absolute()) for i in range(degB)])
        precP = precNewton_mul(A,precA,B,precB)

        #  Using Psi   
        precres = precP*V.newton_polygon()
        quo, precres = Newton_quorem(precres,A0.newton_polygon())
        nwton = sum([ceil(precres(i)) for i in range(d)]) - absA



    # Using Phi ?


    ####################    
    # Lattice estimation
    ####################    
    dA = [x^i for i in range(d)]
    dB = [x^i for i in range(degB+1)]
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
            #scalar = scalar.lift_to_precision(truncprec)
            scalar = scalar.lift()
            Lat[ind] -= scalar * gen



    return abs, nwton, theory, theory2, theory-theory2

def stat_test_factor(NP1,NP2,S,niter,truncprec,nbtest):
    abs=nwton=theory= theory2= theorydiff=0
    for i in range(nbtest):
        abs0, nwton0, theory0, theory20, theorydiff0 = test_factor(NP1,NP2,S,niter,truncprec)
        abs+=abs0
        nwton += nwton0
        theory += theory0
        theory2 += theory20
        theorydiff += theorydiff0
    print("jagged precision")
    print (abs/nbtest).numerical_approx()
    print("Newton precision")
    print (nwton/nbtest).numerical_approx()
    print("lattice precision")
    print (theory/nbtest).numerical_approx(), (theory2/nbtest).numerical_approx(), (theorydiff/nbtest).numerical_approx()
    return (abs/nbtest).numerical_approx(), (nwton/nbtest).numerical_approx(), (theory/nbtest).numerical_approx(), (theory2/nbtest).numerical_approx(), (theorydiff/nbtest).numerical_approx()
        

########################################
# Actual testing
########################################

#####################
#p=2
#####################

p=2
R = Zp(p,prec=100)
prec = 100
SZ.<X> = PolynomialRing(ZZ)
S.<x> = PolynomialRing(R)

#NPtest1 =NewtonPolygon([(0,2),(4,0)])
#NPtest2 =NewtonPolygon([(0,0),(4,3)])

NPtest1 =NewtonPolygon([(0,5),(2,-2),(6,0)])
NPtest2 =NewtonPolygon([(0,0),(4,3),(7,10)])

#stat_test_factor(NPtest1,NPtest2,S,50, 100,1000)
#stat_test_factor(NPtest1,NPtest2,S,50, None,1000)

######

NPtest1 =NewtonPolygon([(0,5),(2,1),(6,0)])
NPtest2 =NewtonPolygon([(0,0),(3,0),(5,10)])

#stat_test_factor(NPtest1,NPtest2,S,50, 100,1000)


#####
NPtest1 =NewtonPolygon([(0,5),(2,1),(3,0)])
NPtest2 =NewtonPolygon([(0,0),(5,-2),(7,10)])

stat_test_factor(NPtest1,NPtest2,S,50, 100,1000)
#stat_test_factor(NPtest1,NPtest2,S,50, None,1000)

########################
#p=11
########################

p=11
R = Zp(p,prec=100)
prec = 100
SZ.<X> = PolynomialRing(ZZ)
S.<x> = PolynomialRing(R)

#NPtest1 =NewtonPolygon([(0,2),(4,0)])
#NPtest2 =NewtonPolygon([(0,0),(4,3)])

NPtest1 =NewtonPolygon([(0,5),(2,-2),(6,0)])
NPtest2 =NewtonPolygon([(0,0),(4,3),(7,10)])

#stat_test_factor(NPtest1,NPtest2,S,50, 100,1000)
#stat_test_factor(NPtest1,NPtest2,S,50, None,1000)

######

NPtest1 =NewtonPolygon([(0,5),(2,1),(6,0)])
NPtest2 =NewtonPolygon([(0,0),(3,0),(5,10)])

#stat_test_factor(NPtest1,NPtest2,S,50, 100,1000)


#####
NPtest1 =NewtonPolygon([(0,5),(2,1),(3,0)])
NPtest2 =NewtonPolygon([(0,0),(5,-2),(7,10)])

##stat_test_factor(NPtest1,NPtest2,S,50, 100,100)
#stat_test_factor(NPtest1,NPtest2,S,50, None,1000)



