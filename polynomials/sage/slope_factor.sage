load("mulmod.sage")
p = 11
R = Zp(11,prec=100)
prec = 100
SZ.<X> = PolynomialRing(ZZ)
S.<x> = PolynomialRing(R)



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
        #my_rem(P,Ai).newton_polygon().plot().show()
        #Ai.newton_polygon().plot().show()
        #print((1-V*Q).mod(Ai))
        #((1-V*Q).mod(Ai).truncate(d)).newton_polygon().plot().show()
    #print(Ai)
    #print(V)
    #Ai.newton_polygon().plot().show()
    #V.newton_polygon().plot().show()
    #P.mod(Ai).truncate(d).newton_polygon().plot().show()
    return Ai

########################################
# Test of slope factorization
########################################

def random_poly_above_NP(NP,R):
    """
	create a polynomial with random coefficients 
        above the Newton polygon     NP.
        highest coefficient is of the form p^l
    """
    p=R.residue_characteristic()
    vert = NP.vertices()
    d = vert[len(vert)-1][0]
    P = p^(ceil(NP(0)))
    for i in range(1,d):
        P += R.random_element()*p^(ceil(NP(i)))*x^i
    P+=p^(ceil(NP(d)))*x^d
    return P
	
def test_factor(NP1,NP2,S,niter):
    """
    create two random polynomials A and B with Newton polygons above NP1 and NP2
    Compute the product AB and factor it by slope factorisation to obtain A
    """
    R=S.base_ring()
    A = random_poly_above_NP(NP1,R)
    B = random_poly_above_NP(NP2,R)
    P = A*B
    vert = NP1.vertices()
    d = vert[len(vert)-1][0]
    A2 = slope_factor(P,d,niter,S)
    return A,B,A2

NPtest1 =NewtonPolygon([(0,2),(4,0)])
NPtest2 =NewtonPolygon([(0,0),(4,3)])


A,B,A2=test_factor(NPtest1,NPtest2,S,50)
P=A*B
#nana=(A-A2).newton_polygon().plot()
Q,R=my_quo_rem(P,A)
Q2,R2=my_quo_rem(P,A2)