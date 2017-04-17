"""
Python implementation of p-adic log (quasi-linear algorithm)
"""

def padiclog(x):
    a = ZZ(x)
    parent = x.parent()
    p = parent.prime()

    if a % p != 1:
        raise NotImplementedError("The argument must be congruent to 1 modulo p")

    prec = x.precision_absolute()
    modulo = p**prec

    N = prec
    while True:
        N2 = prec + floor(log(N)/log(p))
        if N2 == N: break
        N = N2

    trunc = 2; trunc_mod = p*p
    ans = 0

    while True:
        # We extract a small factor f of 1/a
        f = 2 - (a % trunc_mod)
        # We update a
        a *= f

        # We compute -log(f) and add its contribution
        num = [ 1 for i in range(N) ]
        denom = [ i+1 for i in range(N) ]
        step = 1
        h = hpow = 1 - f
        while step < N:
            i = 0
            while i + step < N:
                num[i] = num[i]*denom[i+step] + hpow*num[i+step]*denom[i]
                denom[i] *= denom[i+step]
                i += (step << 1)
            step <<= 1
            hpow = hpow * hpow
        d = num[0].gcd(denom[0])
        num, denom = num[0] // d, denom[0] // d
        _, inv, _ = denom.xgcd(modulo)
        ans += (h * num * inv) % modulo

        # Stop criterium
        if trunc > prec: break

        # We update the truncation level
        trunc_mod = trunc_mod * trunc_mod
        trunc <<= 1
        N >>= 1

    return parent(ans, prec)
