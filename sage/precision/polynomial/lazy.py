from sage.rings.infinity import Infinity

from precision.lazy import lazy_method
from precision.lazy import LazyApproximation


class LazyApproximation_polynomial(LazyApproximation):
    def __init__(self, parent, function, **kwargs):
        if isinstance(function, list):
            from printing import repr_polynomial
            self.repr = lambda: repr_polynomial(function, parent.variable_name())
            kwargs['parenthesis_level'] = 0
        LazyApproximation.__init__(self, parent, function, **kwargs)
        if self._length >= 0:
            self._expected_degree = self._length - 1
        else:
            self._expected_degree = Infinity

    def _add_(self,other,*args,**kwargs):
        res = LazyApproximation._add_(self,other,*args,**kwargs)
        res._expected_degree = max(self._expected_degree, other._expected_degree)
        return res

    def _sub_(self,other,*args,**kwargs):
        res = LazyApproximation._sub_(self,other,*args,**kwargs)
        res._expected_degree = max(self._expected_degree, other._expected_degree)
        return res

    def __neg__(self,*args,**kwargs):
        res = LazyApproximation.__neg__(self,*args,**kwargs)
        res._expected_degree = self._expected_degree
        return res

    def _rmul_(self,*args,**kwargs):
        res = LazyApproximation._rmul_(self,*args,**kwargs)
        res._expected_degree = self._expected_degree
        return res

    def _lmul_(self,*args,**kwargs):
        res = LazyApproximation._lmul_(self,*args,**kwargs)
        res._expected_degree = self._expected_degree
        return res

    def _mul_(self,other,*args,**kwargs):
        res = LazyApproximation._mul_(self,other,*args,**kwargs)
        res._expected_degree = self._expected_degree + other._expected_degree
        return res

    def _pow_(self,exp,*args,**kwargs):
        res = LazyApproximation._pow_(self,exp,*args,**kwargs)
        res._expected_degree = self._expected_degree * exp
        return res
