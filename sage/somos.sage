def somos_naive(a,b,c,d,n):
    if n == 0: return a
    if n == 1: return b
    if n == 2: return c
    for i in range(n-3):
        a, b, c, d = b, c, d, (b*d + c*c) / a
    return d

def change_precision(x, prec):
    if x.precision_absolute() > prec:
        x = x.add_bigoh(prec)
    else:
        x = x.lift_to_precision(prec)
    return x

def somos_stable(a,b,c,d,n):
    if n == 0: return a
    if n == 1: return b
    if n == 2: return c
    v = a.valuation() + b.valuation() + c.valuation() + d.valuation()
    N = min(a.precision_absolute(),
            b.precision_absolute(),
            c.precision_absolute(),
            d.precision_absolute())
    sumv = 0
    maxv = 0
    for i in range(n-3):
        e = (b*d + c*c) / a
        sumv += e.valuation()
        if maxv < e.valuation(): maxv = e.valuation()
        dv = e.valuation() - a.valuation()
        v += dv
        prec = N + v
        prec2 = prec + a.valuation()
        b = change_precision(b, prec2)
        c = change_precision(c, prec2)
        d = change_precision(d, prec2)
        e = (b*d + c*c) / a
        b = change_precision(b, prec)
        c = change_precision(c, prec)
        d = change_precision(d, prec)
        a, b, c, d = b, c, d, e
    print "maxv =", maxv
    return d.add_bigoh(N), sumv

def test(a,b,c,d,n):
    p = a.parent().prime()
    N = min(a.precision_absolute(),
            b.precision_absolute(),
            c.precision_absolute(),
            d.precision_absolute())
    t = cputime()
    r, sumv = somos_stable(a,b,c,d,n)
    t1 = cputime(t)
    print "Time:", t1
    print " r =", r
    print "prec =", N + sumv
    R = Zp(p, N+sumv)
    aa = R(a).lift_to_precision(N+sumv)
    bb = R(b).lift_to_precision(N+sumv)
    cc = R(c).lift_to_precision(N+sumv)
    dd = R(d).lift_to_precision(N+sumv)
    t = cputime()
    rr = somos_naive(aa,bb,cc,dd,n).add_bigoh(N)
    t2 = cputime(t)
    print "Time:", t2
    print "rr =", rr

p = 3
R = Zp(p,100)
a, b, c, d = R(1,20), R(1,20), R(1,20), R(p-1,20)
