p = 2
prec = 200; d = 4
nb = 1000
repeat = 10

R = Zp(p, prec=2*prec)
MS = MatrixSpace(R,d)

def prec_matrix(M):
    return min([ x.precision_absolute() for x in M.list() ])

Basis = [ ]
for i in range(d):
    for j in range(d):
        M = MS(0)
        M[i,j] = 1
        Basis.append(M)

theory = theory2 = practice = 0
for _ in range(repeat):
    Lattice = [ ]
    for i in range(d):
        for j in range(d):
            M = MS(0)
            M[i,j] = 1
            Lattice.append(M)

    M = MS(1)
    for c in range(nb):
        F = MS([ R.random_element().add_bigoh(prec) for _ in range(d^2) ])
        # dM * F + M * dF
        Lat = [ L*F for L in Lattice ] + [ M*L for L in Basis ]
        M = M * F
        #precM = prec_matrix(M)
        #M = MS([ x.add_bigoh(precM) for x in M.list() ])
        Lattice = [ ]
        for i in range(d):
            for j in range(d):
                val = Infinity; index = -1
                for ind in range(len(Lat)):
                    v = Lat[ind][i,j].valuation()
                    if v < val:
                        val = v; index = ind
                if index == -1: raise RuntimeError
                if c == nb-1: theory2 += val
                gen = Lat[index]
                Lattice.append(gen)
                del Lat[index]
                for ind in range(len(Lat)):
                    scalar = Lat[ind][i,j] // gen[i,j]
                    scalar = scalar.lift_to_precision(prec)
                    Lat[ind][i,j] -= scalar * gen[i,j]

    #practice += min([ x.precision_absolute() for x in M.list() ]) - prec
    practice += M[0,0].precision_absolute() - prec
    theory += Lattice[0][0,0].valuation()

print "Theory:", RR(theory/repeat)
print "Theory 2:", RR(theory2/d^2/repeat)
print "Practice:", RR(practice/repeat)
