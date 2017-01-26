def number_diffused(M):
    H = copy(M)
    d = H.ncols()
    n = H.nrows()
    if n < d:
        return Infinity
    number = 0
    for j in range(d):
        val = Infinity; piv = -1
        for i in range(j,n):
            if H[i,j] == 0: continue
            v = H[i,j].valuation()
            if v < val:
                val = v; piv = i
        if val is Infinity: return Infinity
        H.swap_rows(piv,j)
        H.rescale_row(j, ~(H[j,j].unit_part()))
        for i in range(j+1,n):
            H.add_multiple_of_row(i,j, -H[i,j] >> val)

        nb = 0
        for i in range(j):
            v = val - H[i,j].valuation()
            if v > nb: nb = v
        number += nb

    return number


def optimal_jagged(M, diffused=False):
    d = M.nrows()
    if d != M.ncols():
        raise ValueError

    K = M.base_ring()
    S.<X> = PolynomialRing(K)
    C = (X-M).adjoint()

    opt = d * [ Infinity ]

    rows = [ ]
    for i in range(d):
        for j in range(d):
            v = M[j,i].precision_absolute()
            for k in range(d):
                if C[i,j][k] == 0: continue
                w = v + C[i,j][k].valuation()
                if w < opt[k]: opt[k] = w
            if diffused:
                rows.append([ C[i,j][k] << v for k in range(d) ])

    number = 0
    if diffused:
        number = number_diffused(matrix(rows))

    return opt, number


def zealous_jagged(M):
    P = M.charpoly()
    return [ c.precision_absolute() for c in P.list() ]


def stat(d,p=2,prec=200,repeat=10,diffused=False,filename=None):
    if filename is None:
        filename = "results-%s.txt" % d
    Rz = Qp(p,prec=prec)
    Rf = QpFP(p,prec=prec)
    moy_opt = [ 0 ] * d
    moy_zea = [ 0 ] * d
    moy_flo = [ 0 ] * d
    dev_opt = [ 0 ] * d
    dev_zea = [ 0 ] * d
    dev_flo = [ 0 ] * d
    moy_diff = 0; dev_diff = 0
    count = 0
    import io
    fh = io.open(filename,'a')
    fh.write(u"	Opt.	CR	FP\n")
    for _ in range(repeat):
        Mz = random_matrix(Rz,d,d)

        optimal, number = optimal_jagged(Mz, diffused=diffused)
        if Infinity in optimal or number is Infinity:
            continue

        Pz = Mz.charpoly(algorithm="df")
        zealous = [ Pz[i].precision_absolute() for i in range(d) ]

        diffprec = max([ optimal[i] - zealous[i] for i in range(d) ])
        Rz2 = Qp(p, prec=prec+diffprec)
        Mz2 = matrix(d,d,[ Rz2(c).lift_to_precision() for c in Mz.list() ])
        Pz2 = Mz2.charpoly()

        Mf = Mz.change_ring(Rf)
        Pf = Mf.charpoly(algorithm="df")
        floating = [ min(optimal[i], (Pz2[i]-Rz2(Pf[i])).valuation()) for i in range(d) ]

        for i in range(d):
            v = Pz2[i].valuation()
            opt = prec - (optimal[i] - v)
            zea = prec - (zealous[i] - v)
            flo = prec - (floating[i] - v)
            fh.write(u"X^%s	%s	%s	%s\n" % (i, opt, zea, flo))
            moy_opt[i] += opt
            moy_zea[i] += zea
            moy_flo[i] += flo
            dev_opt[i] += opt**2
            dev_zea[i] += zea**2
            dev_flo[i] += flo**2
        if diffused:
            moy_diff += number
            dev_diff += number**2
            fh.write(u"Diff.	%s\n" % number)
        fh.write(u"\n")
        fh.flush()

        count += 1
        if count % 10 == 0: print count

    fh.write(u"Average\n")
    for i in range(d):
       moy_opt[i] = RR(moy_opt[i]/count)
       moy_zea[i] = RR(moy_zea[i]/count)
       moy_flo[i] = RR(moy_flo[i]/count)
       fh.write(u"X^%s	%s	%s	%s\n" % (i,moy_opt[i], moy_zea[i], moy_flo[i]))
    if diffused:
       moy_diff = RR(moy_diff/count)
       fh.write(u"Diff:	%s\n" % moy_diff)

    fh.write(u"\nDeviation\n")
    for i in range(d):
       dev_opt[i] = sqrt(RR(dev_opt[i]/count) - moy_opt[i]**2)
       dev_zea[i] = sqrt(RR(dev_zea[i]/count) - moy_zea[i]**2)
       dev_flo[i] = sqrt(RR(dev_flo[i]/count) - moy_flo[i]**2)
       fh.write(u"X^%s	%s	%s	%s\n" % (i,dev_opt[i], dev_zea[i], dev_flo[i]))
    if diffused:
       dev_diff = sqrt(RR(dev_diff/count) - moy_diff**2)
       fh.write(u"Diff:	%s\n" % dev_diff)

    fh.write(u"\nRatio of success: %s\n" % RR(count/repeat))

    fh.write(u"\n----\n\n")
    fh.close()
