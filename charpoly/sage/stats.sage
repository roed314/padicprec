def optimal_jagged(M):
    d = M.nrows()
    if d != M.ncols():
        raise ValueError

    K = M.base_ring()
    S.<X> = PolynomialRing(K)
    C = (X-M).adjoint()

    opt = d * [ Infinity ]

    for i in range(d):
        for j in range(d):
            v = M[j,i].precision_absolute()
            for k in range(d):
                if C[i,j][k] == 0: continue
                w = v + C[i,j][k].valuation()
                if w < opt[k]: opt[k] = w

    return opt


def zealous_jagged(M):
    P = M.charpoly()
    return [ c.precision_absolute() for c in P.list() ]


def test(d,p=2,prec=200,repeat=1,filename=None):
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
    count = 0
    import io
    fh = io.open(filename,'a')
    fh.write(u"	Opt.	CR	FP\n")
    for _ in range(repeat):
        Mz = random_matrix(Rz,d,d)

        optimal = optimal_jagged(Mz)
        if Infinity in optimal:
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
            fh.write(u"X^%s	%s	%s	%s\n" % (i,optimal[i]-v, zealous[i]-v, floating[i]-v))
            moy_opt[i] += optimal[i] - v
            moy_zea[i] += zealous[i] - v
            moy_flo[i] += floating[i] - v
            dev_opt[i] += (optimal[i] - v)**2
            dev_zea[i] += (zealous[i] - v)**2
            dev_flo[i] += (floating[i] - v)**2
        fh.write(u"\n")
        fh.flush()

        count += 1

    fh.write(u"Average\n")
    for i in range(d):
       moy_opt[i] = RR(moy_opt[i]/count)
       moy_zea[i] = RR(moy_zea[i]/count)
       moy_flo[i] = RR(moy_flo[i]/count)
       fh.write(u"X^%s	%s	%s	%s\n" % (i,moy_opt[i], moy_zea[i], moy_flo[i]))

    fh.write(u"Deviation\n")
    for i in range(d):
       dev_opt[i] = sqrt(RR(dev_opt[i]/count) - moy_opt[i]**2)
       dev_zea[i] = sqrt(RR(dev_zea[i]/count) - moy_zea[i]**2)
       dev_flo[i] = sqrt(RR(dev_flo[i]/count) - moy_flo[i]**2)
       fh.write(u"X^%s	%s	%s	%s\n" % (i,dev_opt[i], dev_zea[i], dev_flo[i]))

    fh.write(u"Ratio of success: %s\n" % RR(count/repeat))

    fh.write(u"\n----\n\n")
    fh.close()

    return moy_opt, dev_opt, moy_zea, dev_zea, moy_flo, dev_flo, RR(count/repeat)
