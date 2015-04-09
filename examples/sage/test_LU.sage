#Computation of low(A) and up(A)#

def low(M):
    n = M.nrows()
    m = M.ncols()
    A = copy(M)
    for i in range(n):
        for j in range(i,m):
            A[i,j] = 0
    return A

def up(M):
    n = M.nrows()
    m = M.ncols()
    A = copy(M)
    for i in range(n):
        for j in range(i):
            A[i,j] = 0
    return A

def printval(M,p):
    n = M.nrows()
    m = M.ncols()
    Mval = MatrixSpace(QQ,n,m)(0)
    for i in range(n):
        for j in range(m):
            if M[i,j].valuation(p) ==+Infinity:
                Mval[i,j] = -42
            else:
                Mval[i,j] =M[i,j].valuation(p)
    print(Mval)
    return Mval


def test_DM_LU(d,p,prec):

    R = Zp(p, prec=2*prec)
    MS = MatrixSpace(R,d)

    MEX = MatrixSpace(QQ,d)

    #Basis for MEX

    Basis = [ ]
    listcouple = [ ]
    for ii in range(d):
        for jj in range(d):
            Mb = MEX(0)
            Mb[ii,jj] = 1
            Basis.append(Mb)
            listcouple.append([ii,jj])
    


    # Construction d'une matrice aleatoire #

    Mzp = MS.random_element()
    M = MEX(0)
    for i in range(d):
        for j in range(d):
            M[i,j] = (Mzp[i,j]).lift()
            
    #print("################")
    #print("forme de M")
    #printval(M,p)



    lu = M.LU(pivot ='nonzero')
    L = lu[1]
    U = lu[2]
    
    #print("################")
    #print("forme de L et U")
    #printval(L,p)
    #printval(U,p)
    #print("forme de L^-1 et U^-1")
    #printval(L^(-1),p)
    #printval(U^(-1),p)

    #Calcul de dM#

    MEX2 = MatrixSpace(QQ,d^2)

    DM = MEX2(0)



    for k in range(d^2):
        u = listcouple[k][0]
        v = listcouple[k][1]
        dM = Basis[k]
        meuv =L*low(L^(-1)*dM*U^(-1))+up(L^(-1)*dM*U^(-1))*U
        for kk in range(d^2):
            uu = listcouple[kk][0]
            vv = listcouple[kk][1]
            DM[kk,k] = meuv[uu,vv]
    
    print("####################")
    print("forme de la differentielle")
    printval(DM,p)

            
    #test de la differentielle
    perturb = p^prec*MS.random_element()
    Mzp2=Mzp+perturb
    M2 = MEX(0)
    dMth = MEX(0)
    for iii in range(d):
        for jjj in range(d):
            M2[iii,jjj] = (Mzp2[iii,jjj]).lift()
            dMth = (perturb[iii,jjj]).lift()
    
    #print("####################")
    #print("les modifies")
    #printval(M-M2,p)
    #printval(M2,p)
    
    DLth = L*low(L^(-1)*dMth*U^(-1))
    DUth = up(L^(-1)*dMth*U^(-1))*U
    
    lu2 = M2.LU(pivot ='nonzero')
    L2 = lu2[1]-L
    U2 = lu2[2]-U
    
    #print("####################")
    #print("Diff effective de L")
    #printval(L2,p)
    #print("Diff effective de U")
    #printval(U2,p)

    #print("####################")
    #print("Diff theorique de L")
    #printval(DLth,p)
    #print("Diff theorique de U")
    #printval(DUth,p)
    
    Ldiff = L2-DLth
    Udiff = U2-DUth
     
    
    print("################")
    print("test de dL")
    printval(Ldiff,p)
    print("################")
    print("test de dU")
    printval(Udiff,p)    



    #Calcul de la valuation pour la projection#
    valproj = 0
    for l in range(d^2):
        valproj += min([DM[l,ll].valuation(p) for ll in range(d^2)])

    #Calcul de la perte de precision pour les lattices#
    vallatt = DM.det().valuation(p)




test_DM_LU(3,2,75)