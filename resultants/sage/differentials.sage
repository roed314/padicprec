p = 3
prec = 100
degree = 15

K = Qp(p,prec)
R = K.integer_ring()
S.<x> = PolynomialRing(K)

P = S( [ R.random_element().add_bigoh(prec) for _ in range(degree+1) ] )
Q = S( [ R.random_element().add_bigoh(prec) for _ in range(degree+1) ] )
_,U,V = P.xgcd(Q)

rows = [ ]

for i in range(degree):
    # dA = x^i, dB = 0
    dU = x^i*U^2 % Q
    dV = x^i*U*V % P
    row = [ dU[j] for j in range(degree) ] + [ dV[j] for j in range(degree-1) ]
    rows.append(row)

    # dA = 0, dB = x^i
    dU = x^i*U*V % Q
    dV = x^i*V^2 % P
    row = [ dU[j] for j in range(degree) ] + [ dV[j] for j in range(degree-1) ]
    rows.append(row)

M = matrix(rows)
MM = M.__copy__()

for i in range(2*degree-1):
    val = 1000
    for ii in range(i,2*degree):
        for jj in range(i,2*degree-1):
            v = MM[ii,jj].valuation()
            if v < val:
                val = v
                pivot = (ii,jj)
    print val
    (ii,jj) = pivot
    MM.swap_rows(ii,i)
    MM.swap_columns(jj,i)
    inv = ~MM[i,i]
    for ii in range(i+1,2*degree):
        MM.add_multiple_of_row(ii, i, -inv*MM[ii,i])
    for jj in range(i+1,2*degree-1):
        MM.add_multiple_of_column(jj, i, -inv*MM[i,jj])

print "--";

loss = 0
for i in range(degree):
    d = min(U[i].precision_absolute(), V[i].precision_absolute()) - prec
    if d < loss: loss = d
print loss
