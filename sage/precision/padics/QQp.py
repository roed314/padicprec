from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.integer import Integer
from sage.rings.infinity import Infinity

from sage.rings.ring import CommutativeRing
from sage.rings.morphism import RingHomomorphism_coercion

from precision.parent_inexact import Parent_inexact
from approximation import QQpApprox
from precision.bigoh import ValuationBigOh


class pAdicCoersion(RingHomomorphism_coercion):
    def __init__(self,A,B):
        RingHomomorphism_coercion.__init__(self, A.Hom(B), check=False)

    def _call_(self,x):
        from QQp_element import QQp_element
        return QQp_element(self.codomain(), x)

    def _call_with_args(self,x,args,kwargs):
        from QQp_element import QQp_element
        return QQp_element(self.codomain(), x, *args, **kwargs)


class QQp(Parent_inexact,CommutativeRing):  # Should create RingInexact
    def __init__(self, p):
        if not p.is_prime():
            raise ValueError("p must be a prime number")
        self._p = p
        Parent_inexact.__init__(self, self, QQpApprox(p))
        from ZZp import ZZp
        injectQQ = pAdicCoersion(QQ, self)
        injectZp = pAdicCoersion(ZZp(p), self)
        self._populate_coercion_lists_(coerce_list=[injectQQ,injectZp])

    def _element_constructor_(self, x, *args, **kwargs):
        from QQp_element import QQp_element
        if isinstance(x,int): x = ZZ(x)
        return QQp_element(self, x, *args, **kwargs)

    def p(self):
        return self._p

    def dimension(self):
        return 1

    def __repr__(self):
        return "%s-adic Field" % self._p

    def is_commutative(self):
        return True

    def uniformizer_name(self):
        return str(self._p)

    def uniformizer(self,prec=Infinity):
        from QQp_element import QQp_element
        return QQp_element(self, self._p, prec)
