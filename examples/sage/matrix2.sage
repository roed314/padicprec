p = 2
prec = 1000; d = 3
nb = 500
repeat = 100

R = Zp(p, prec=2*prec)
MS = MatrixSpace(R,d)

def prec_matrix(M):
    return min([ x.precision_absolute() for x in M.list() ])

results = [ ]
for _ in range(repeat):
    M = MS(1)
    for c in range(nb):
        F = MS([ R.random_element().add_bigoh(prec) for _ in range(d^2) ])
        M = M * F
    vp1 = min([ x.valuation() for x in M.list() ])
    vp2 = M.determinant().valuation()/d
    results.append([vp1,vp2])
