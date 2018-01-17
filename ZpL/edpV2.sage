##################
#      4.1       #
##################

def edp_simple(g,h,N):
    xe=g.parent().gen()
    u =xe
    for i in range(N):
        print i
        uprime = u.derivative()
        temp=h(u)*((uprime/h(u)-g).integral())
        u = (u -temp)
        u=sum([u[k]*xe^k for k in range(2^(i+1))])
        print u 
        print "#################"
    return u

pprec=250
xprec = 2^5
A=QpLC(2,pprec,print_mode='terse')
B=Zp(2,pprec)
R.<x> = PolynomialRing(A)
St.<t> = PowerSeriesRing(A,xprec)
Stint.<ti>=PowerSeriesRing(B,xprec)
h0=1+t+t^3
y=t+t^2*St(Stint.random_element())
#y=t+t^2-2*t^3+3*t^4
g0= St(St(y.derivative())/St(h0(y)))




S=edp_simple(g0,h0,5)
N=2^5
(S.derivative()).truncate_powerseries(N)-(g0*h0(S)).truncate_powerseries(N)
[(y-S)[k].valuation() for k in range(N)] 
sage: [(y-S)[k].valuation() for k in range(N)] 


A2=Qp(2,pprec,print_mode='terse')
St2.<t2> = PowerSeriesRing(A2,xprec)
h02=sum([h0[k].lift()*t2^k for k in range(xprec)])
y2=sum([y[k].lift()*t2^k for k in range(xprec)])
#y=t+t^2-2*t^3+3*t^4
g02= sum([g0[k].lift()*t2^k for k in range(xprec)])
S2=edp_simple(g02,h02,5)
N=2^5
[(y2[k].lift()-S2[k].lift()).valuation(2) for k in range(N)] 



##################
#      5         #
##################







    
def Newton_DE_quad(R,g,h,N,m):
    alpha = 0
    beta = 1
    d = 2
    U = 1/beta
    J = 1
    V = 1
    #S = (alpha + beta*x+x^2*(g.derivative()+h.derivative()*beta^3)/(4*beta)).truncate(1)
    S = R.gen()
    for j in range(1,N):
        print j
        Sprime = S.derivative()
        compo= h(S)
        U = (U*(2-Sprime*U)).truncate(2^j) 
        V = R((V+J*compo*(2-V*J))/2).truncate(2^j) 
        J = (J*(2-V*J)).truncate(2^j) 
        S = R(S +V*(((U*J/2)*(g*compo-Sprime^2)).integral())).truncate(2^(1+j)) 
        d = 2*d
        print(S)
    return S

#QpL
A = QpLC(3,50,print_mode='terse')
prec = A.precision()
#A = QQ
N=4
R.<x> = PolynomialRing(A)

m=67

g=1+(1/4)*x^2+x^6
h=1+((m^2)/4)*x^2+m^6* x^6 
S= Newton_DE_quad(R,g,h,6,m) 




#Qp
A = Qp(3,50,print_mode='terse')
#A = QQ
N=4
R.<x> = PolynomialRing(A)

m=67

g=1+(1/4)*x^2+x^6
h=1+((m^2)/4)*x^2+m^6* x^6 
S= Newton_DE_quad(R,g,h,6,m) 
