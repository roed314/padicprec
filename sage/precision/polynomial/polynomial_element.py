from sage.rings.infinity import Infinity

from precision.element_inexact import RingElement_inexact
from precision.bigoh import ExactBigOh

from precision.exception import PrecisionError
from precision.limits import limit_lazyprecision


class Polynomial_inexact(RingElement_inexact):
    def _normalize(self):
        pass

    def variable_name(self):
        return self.parent().variable_name()

    def valuation_bigoh(self,model=None,lazylimit=Infinity):
        from bigoh import FlatIntervalBigOh_polynomial, JaggedBigOh_polynomial, NewtonBigOh_polynomial, LatticeBigOh_polynomial
        deg = self.degree(minimal=False)
        valuations = [ self[i].valuation(minimal=True,lazylimit=lazylimit) for i in range(deg+1) ]
        valuations.append(Infinity)
        parentprec = self._precision.parent()
        if model is None:
            model = self.precision_model()
        if model is FlatIntervalBigOh_polynomial:
            val = min(valuations)
            first = len(valuations)
            for i in range(len(valuations)):
                if valuations[i] < Infinity:
                    first = i
                    break
            last = 0
            for i in range(len(valuations), 0, -1):
                if valuations[i-1] < Infinity:
                    last = i
                    break
            return FlatIntervalBigOh_polynomial(parentprec, val, last, first)
        elif model is JaggedBigOh_polynomial:
            return JaggedBigOh_polynomial(parentprec, valuations)
        elif model is NewtonBigOh_polynomial or model is ExactBigOh:
            return NewtonBigOh_polynomial(parentprec, valuations)
        else:
            raise TypeError("invalid model")

    def degree(self,minimal=True,lazylimit=Infinity):
        deg = -1
        for d in range(self._approximation._expected_degree,-1,-1):
            if not self[d].is_zero(lazylimit=lazylimit): 
                deg = d
                break
        if minimal:
            return deg
        else:
            return max(deg, self.precision().last()-1)

    def add_bigoh(self,prec,degree=None):
        from bigoh import FlatIntervalBigOh_polynomial
        parentprec = self.parent().precision()
        precision = self._precision
        try:
            precscal = parentprec._base_precision(prec)
            if degree == None:
                degree = self.degree(minimal=False)
            precision += FlatIntervalBigOh_polynomial(parentprec,precscal,degree+1)
        except PrecisionError:
            precision += parentprec(prec)
        return self._new_fast(self.parent(), self._approximation, precision)

    def list(self):
        return [ self[i] for i in range(self.degree(minimal=False,lazylimit=limit_lazyprecision)+1) ]

    def _repr_(self):
        approximation = self._approximation
        length = max(approximation.length(), self.precision().last())
        if length >= 0:
            coeffs = [ self[i] for i in range(length) ]
            from printing import repr_polynomial
            return repr_polynomial(coeffs, self.variable_name())
        else:
            return RingElement_inexact._repr_(self)

    def newton_polygon(self,minimal=False):
        deg = self.degree(not minimal)
        vertices = [ (i,self[i].valuation(minimal)) for i in range(deg+1) ]
        from sage.geometry.newton_polygon import NewtonPolygon
        return NewtonPolygon(vertices)

    def newton_slopes(self, minimal=False, repetition=True):
        slopes = [ -s for s in self.newton_polygon(minimal=minimal).slopes(repetition=repetition) ]
        slopes.reverse()
        return slopes
