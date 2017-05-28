# log(a(x)) = x + x^2/2
# 
#     a_0 = 1 + x
# a_(i+1) = a_i * (1 + x + x^2/2 - log(a_i))
# 
# Posons l_i = log(a_i)
#          y = x + x^2/2
#        b_i = y - l_i
# 
# a_0 = 1+x
# l_0 = log(1+x)
# b_0 = x + x^2/2 - log(1+x)
# 
# a_(i+1) = a_i * (1 + y - l_i)
# l_(i+1) = l_i + log(1 + y - l_i)
#         = l_i + (y - l_i) - (y - l_i)^2/2 + ...
#         = y - (y - l_i)^2/2 + ...
# 
# a_(i+1) = a_i * (1 + b_i)
# b_(i+1) = log(1 + b_i) - b_i
#
#
# Version tronquée 
#
# a_0 = 1+x
# l_0 = log(1+x)
# b_0 = 0
# 
# a_(i+1) = a_i * (1 + b_i)
# l_(i+1) = l_i + log(1 + b_i)
# b_(i+1) = y - l_(i+1)   mod  p^(...)


def artin_hasse(x, ans=None):
    R = x.parent()
    prec = x.precision_absolute()
    p = R.prime()

    y = xp = x; i = 1
    while True:
        xp = xp**p
        if xp.valuation() >= prec + i: break
        y += (xp >> i)
        i += 1

    a = (1+x).add_bigoh(2).lift_to_precision()
    l = log(a)
    b = (y-l).add_bigoh(3).lift_to_precision()
    trunc = 2
    while trunc < prec:
        trunc = 2*trunc - 1
        a *= 1+b
        l += log(1+b)
        b = (y-l).add_bigoh(2*trunc-1).lift_to_precision()
        if ans is not None:
            print trunc, (ans-a).valuation()
    return a
