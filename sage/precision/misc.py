from sage.rings.infinity import Infinity

def ceil_infty(x):
    if x is Infinity:
        return x
    else:
        return x.ceil()
