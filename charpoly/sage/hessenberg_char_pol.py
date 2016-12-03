def Hessenberg(Mat):
    r"""
    
    Computes the Hessenberg form of the matrix M,
    along with the passage matrix P

    
    INPUT:

    - ``Mac`` -- a matrix
    
    OUTPUT:

    H, a matrix under Hessenberg form, P the passage matrix, Pinv its inverse.

    EXAMPLES::

        sage: M = Matrix(Qp(2,50),4,4,[1,2,3,4,-5,-6,7,8,9,10,-11,12,13,14,15,-16])
        sage: H,P,Pinv = Hessenberg(M)
        sage: H
        [  0   1   2   0]
        [  1   0  -4  -3]
        [  0   0   1 3/4]
        [  0   0   0   1]
        sage: P
        [1, 0, 2, 3]
        sage: P*M*Pinv

    """
    M = copy(Mat)
    n = M.nrows()
    m = M.ncols()
    K = M.base_ring()
    R = K.integer_ring()

    #Beware: P could also be defined at infinite precision, but this is not our choice today
    P = MatrixSpace(K,n,n)(1)
    Pinv = MatrixSpace(K,n,n)(1)

    for j in range(m-1):
        i=j+1
        imax=-1
        #cmax= M[i,j]
        while i<n:
            if M[i,j] <>0:
                if imax==-1:
                    imax=i
                    cmax= M[i,j]
                else:
                    if M[i,j].valuation()<cmax.valuation():
                        imax=i
                        cmax= M[i,j]
            i=i+1
        if imax>=0:
            #M.rescale_row(i, M[i,j]^(-1), start_col=j)
            M.swap_rows(j+1,imax)
            M.swap_columns(j+1,imax)

            P.swap_rows(j+1,imax)
            Pinv.swap_columns(j+1,imax)

            for l in range(j+2,n):
                coeff = (M[l,j]*M[j+1,j]**(-1)).lift()
                M.add_multiple_of_row(l, j+1, -coeff , start_col=j)
                P.add_multiple_of_row(l, j+1, -coeff)

                M.add_multiple_of_column(j+1, l, coeff)
                Pinv.add_multiple_of_column(j+1, l, coeff)                
 
    return M,P,Pinv

def IdXA(A):
    r"""
    
    Computes the Id-XA for A a square matrix

    
    INPUT:

    - ``A`` -- a square matrix
    
    OUTPUT:

    Id-XA

    EXAMPLES::

        sage: M = Matrix(QQ,4,4,[1,2,3,4,-5,-6,7,8,9,10,-11,12,13,14,15,-16])
        sage: B = IdXA(M)
        sage: B
        [  0   1   2   0]
        [  1   0  -4  -3]
        [  0   0   1 3/4]
        [  0   0   0   1]


    """
    n = A.nrows()
    K = A.base_ring()
    P = PowerSeriesRing(K, name='X', default_prec = n+1)
    X = P.gen()
    MatXspace = MatrixSpace(P,n,n)
    B = MatXspace(1)-X*MatXspace(A)
    return B


def resolvent_and_char_pol_of_Hessenberg(A):
    r"""
    
    Computes the comatrix of Id-XA for A a square matrix of dimension n,
    from the resolvent of A, and the reciprocal of its characteristic polynomial

    
    INPUT:

    - ``A`` -- a square matrix
    
    OUTPUT:

    Com(Id-XA), det(Id-XA)

    EXAMPLES::

        sage: A = Matrix(QQ,4,4,[1,2,3,4,-5,-6,7,8,0,10,-11,12,0,0,15,-16])
        sage: B = IdXA(A)
        sage: C = resolvent_and_char_pol_of_Hessenberg(A)
        [  0   1   2   0]
        [  1   0  -4  -3]
        [  0   0   1 3/4]
        [  0   0   0   1]
        sage: B*C

    """
    n = A.nrows()
    M = IdXA(A)
    P = M.base_ring()
    x = P.gen()
    Minv = M.parent()(1)
    Qinv = M.parent()(1)
    Pinv = M.parent()(1)
    Q = M.parent()(1)
    P = M.parent()(1)


    #column elimination
    
    for i in range(n):
        #ci_inv = M[i,i].inverse_mod(x^(n+1))
        ci_inv = M[i,i]**(-1)
        for j in range(i+1,n):    
            coeff = M[i,j]*ci_inv
            # Beware, we did not put any start_row or end row...
            M.add_multiple_of_column(j,i , -coeff)
            Q.add_multiple_of_column(j, i, -coeff)
            Qinv.add_multiple_of_row(i,j,coeff)



    
    #row elimination


    ci_inv = []

    for i in range(n-1):
        #ci_inv.append(M[i,i].inverse_mod(x^(n+1)))
        ci_inv.append(M[i,i]**(-1))
        coeff = M[i+1,i]*ci_inv[i]
        M.add_multiple_of_row(i+1,i , -coeff)
        P.add_multiple_of_row(i+1, i, -coeff)
        Pinv.add_multiple_of_column(i,i+1,coeff)

    #ci_inv.append(M[n-1,n-1].inverse_mod(x^(n+1)))
    ci_inv.append(M[n-1,n-1]**(-1))




    # reading the reciprocal char pol
    # Beware: precision is not done properly here

    rec_chi = prod([M[i,i] for i in range(n)])
    R = PolynomialRing(M.parent().base_ring().base_ring(), name='X')
    chi = R([ rec_chi[n-i] for i in range(n+1) ])


    #final normalization

    for i in range(n):
        Pinv.rescale_row(i,M[i,i])
        P.rescale_row(i,ci_inv[i])


    #final computation for comatrix
    Minv = Q*P
    comat = Minv*rec_chi



    return comat, chi
        
          

#M = Matrix(Qp(2,50),4,4,[1,2,3,4,-5,-6,7,8,9,10,-11,12,13,14,15,-16])
#H,P,Pinv = Hessenberg(M)
#M = Matrix(QQ,4,4,[1,2,3,4,-5,-6,7,8,9,10,-11,12,13,14,15,-16])
#b = IdXA(M)
#
#A = Matrix(QQ,4,4,[1,2,3,4,-5,-6,7,8,0,10,-11,12,0,0,15,-16])
#B = IdXA(A)
#C,chi = resolvent_and_char_pol_of_Hessenberg(A)
#print(B*C)
