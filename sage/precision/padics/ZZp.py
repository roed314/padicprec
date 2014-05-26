from sage.rings.integer_ring import ZZ
from sage.rings.integer import Integer
from sage.rings.infinity import Infinity
from sage.rings.ring import Ring

from sage.rings.morphism import RingHomomorphism_coercion

from precision.parent_inexact import Parent_inexact
from precision.bigoh import ValuationBigOh
from approximation import QQpApprox


class pAdicCoersion(RingHomomorphism_coercion):
    def __init__(self,A,B):
        RingHomomorphism_coercion.__init__(self, A.Hom(B), check=False)

    def _call_(self,x):
        from ZZp_element import ZZp_element
        return ZZp_element(self.codomain(), x)

    def _call_with_args(self,x,args,kwargs):
        from ZZp_element import ZZp_element
        return ZZp_element(self.codomain(), x, *args, **kwargs)


class ZZp(Parent_inexact,Ring):  # Should create RingInexact
    def __init__(self, p, capped_relative=None, capped_absolute=None, truncate_precision=20):
        if not p.is_prime():
            raise ValueError("p must be a prime number")
        self._p = p
        Parent_inexact.__init__(self, self, QQpApprox(p), capped_relative, capped_absolute)
        self._truncate_precision = truncate_precision
        injectZZ = pAdicCoersion(ZZ, self)
        self._populate_coercion_lists_(coerce_list=[injectZZ])
        from lazy import LazyApproximation_padics
        self._lazy_class = LazyApproximation_padics

    def _element_constructor_(self, x, *args, **kwargs):
        from ZZp_element import ZZp_element
        if isinstance(x,int): x = ZZ(x)
        return ZZp_element(self, x, *args, **kwargs)

    def dimension(self):
        return 1

    def __repr__(self):
        s = "%s-adic Ring" % self._p
        capped_relative = self.capped_relative()
        if capped_relative is Infinity:
            word = "with"
        else:
            s += " with capped relative precision %s" % capped_relative
            word = "and"
        capped_absolute = self.capped_absolute().flat_below()
        if capped_absolute is not Infinity:
            if capped_absolute > capped_relative:
                s += " %s capped absolute precision %s" % (word, capped_absolute)
            else:
                s += " (%s capped absolute precision %s)" % (word, capped_absolute)
        return s

    def is_commutative(self):
        return True

    def uniformizer_name(self):
        return str(self._p)

    def uniformizer(self,prec=Infinity):
        from ZZp_element import ZZp_element
        return ZZp_element(self, self._p, prec)

    def random_element(self,prec=Infinity,precrandom=None):
        from ZZp_element import ZZp_element
        if precrandom is None:
            if prec is not Infinity:
                precrandom = prec
            else:
                precrandom = self._prec
        modulo = self._p ** precrandom
        approximation = ZZ.random_element(modulo)
        precision = self.precision()(prec)
        return ZZp_element(self, approximation, precision)
