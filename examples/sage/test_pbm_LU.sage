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
    
def test_LU_minimal(d,p,prec):
    R = Zp(p, prec=2*prec)
    MS = MatrixSpace(R,d)

    MEX = MatrixSpace(QQ,d)

    Mzp = MS.random_element()
    M = MEX(0)
    for i in range(d):
        for j in range(d):
            M[i,j] = (Mzp[i,j]).lift()
    
    print("################")
    print("################")
    print("################")
    #print("forme de M")
    #printval(M,p)
    
    lu = M.LU(pivot ='nonzero')
    L = lu[1]
    U = lu[2]
    
    #print("################")
    #print("forme de L")
    #printval(L,p)
    #print("forme de U")
    #printval(U,p)
    
    perturb = p^prec*MS.random_element()
    Mzp2=Mzp+perturb
    M2 = MEX(0)
    dMth = MEX(0)
    for iii in range(d):
        for jjj in range(d):
            M2[iii,jjj] = (Mzp2[iii,jjj]).lift()
            dMth[iii,jjj] = (perturb[iii,jjj]).lift()
    
    #print("####################")
    #print("la perturbation")
    #printval(dMth,p)    
    #print("les modifies")
    #printval(M2-M,p)
    #printval(M2,p)
    
    lu2 = M2.LU(pivot ='nonzero')
    L2 = lu2[1]-L
    U2 = lu2[2]-U
    
    #print("####################")
    #print("Nouveau L")
    #printval(lu2[1],p)
    #print("Nouveau U")
    #printval(lu2[2],p)  
    
    print("####################")
    print("Diff effective de L")
    printval(L2,p)
    print("Diff effective de U")
    printval(U2,p)
    
    return M,M2
    
def repeat_test(nbrepeat,d,p,prec):
    for i in range(nbrepeat):
        test_LU_minimal(d,p,prec)

d=2;p=2;prec=75;
repeat_test(4,d,p,prec)
