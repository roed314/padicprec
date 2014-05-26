from sage.rings.infinity import Infinity
from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
from sage.rings.polynomial.polynomial_element import Polynomial_generic_dense

from precision.approximation import Approximation

class PolynomialRingApprox(PolynomialRing_general):
    def __init__(self, *args, **kwargs):
        PolynomialRing_general.__init__(self, element_class=PolynomialApprox, *args, **kwargs)
        from lazy import LazyApproximation_polynomial
        self._lazy_class = LazyApproximation_polynomial
        self._zero = self(0)
        self._lazy_zero = self._lazy_class(self, self._zero)

    def indices_basis(self):
        pass

    def dimension(self):
        return Infinity


class PolynomialApprox(Polynomial_generic_dense, Approximation):
    def __init__(self, parent, x, check=True, is_gen=False, construct=False):
        Polynomial_generic_dense.__init__(self, parent, x, check=check, is_gen=is_gen, construct=construct)
        Approximation.__init__(self,parent)
        self._expected_degree = self.degree()
        self._length = self._expected_degree + 1

    def truncate(self,workprec):
        return self.__class__(self.parent(), [ c.truncate(workprec) for c in self.list() ])

    def parenthesis_level(self):
        return 0

    def _getitem_by_num(self,i):
        return self[i]

    # Grrr
    def _add_(self,other):
        res = Polynomial_generic_dense._add_(self,other)
        Approximation.__init__(res,res.parent())
        res._expected_degree = res.degree()
        res._length = res._expected_degree + 1
        return res
    def _mul_(self,other):
        res = Polynomial_generic_dense._mul_(self,other)
        Approximation.__init__(res,res.parent())
        res._expected_degree = res.degree()
        res._length = res._expected_degree + 1
        return res
