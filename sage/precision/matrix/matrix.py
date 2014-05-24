from operator import mul

from sage.rings.infinity import Infinity
from sage.structure.element import Element

from precision.element_inexact import Element_inexact, ModuleElement_inexact
from precision.bigoh import BigOh, ExactBigOh
from precision.exception import PrecisionError
from precision.grouping.mode import LazyMode_Disable

from precision.limits import limit_lazyprecision


class Matrix_inexact(ModuleElement_inexact):
    def __init__(self, parent, x, prec=None, workprec=None, group=None):
        if parent.is_square():
            base = parent.base_ring()
            if isinstance(x,Element) and x.parent() is base:
                c = x
            else:
                try:
                    c = base(x)
                except (TypeError,ValueError,PrecisionError):
                    c = None
            if c is not None:
                n = parent.nrows()
                x = (n**2) * [ base(0) ]
                for i in range(0,n**2,n+1):
                    x[i] = c
        ModuleElement_inexact.__init__(self, parent, x, prec, workprec, group)

    def is_square(self):
        return self.parent().is_square()

    def __mul__(self,other):
        parent = self.parent()
        if isinstance(other,Matrix_inexact) and parent.base_ring() is other.parent().base_ring():
            group = self.group()
            if other.group() != group:
                raise TypeError("Can multiply only two matrices in the same group")
            return self._matrix_times_matrix_(other)
        else:
            from sage.structure.element import get_coercion_model
            return get_coercion_model().bin_op(self,other,mul)

    def _matrix_times_matrix_(self,other):
        if self.ncols() != other.nrows():
            raise TypeError("incompatible dimensions")
        from matrix_space import MatrixSpace_inexact
        parent = MatrixSpace_inexact(self.parent().base_ring(), self.nrows(), other.ncols())
        selfval = self.valuation_bigoh(lazylimit=limit_lazyprecision)
        otherval = other.valuation_bigoh(lazylimit=limit_lazyprecision)
        precision = self._precision * otherval + selfval * other._precision
        if self._group.lazy_mode() is not LazyMode_Disable:
            selfv = selfval.flat_below()
            otherv = otherval.flat_below()
            approximation = self._approximation._mul_(other._approximation, workprec_shifts=[-otherv,-selfv], expected_valuation = selfv+otherv)
        else:
            approximation = self._approximation._mul_(other._approximation)
        return self._new_fast(parent, approximation, precision)

    def __pow__(self,exp,ignored=None):
        if not self.is_square():
            raise TypeError("The matrix must be square")

        # Copied from precision.element_inexact.RingElement_inexact
        parent = self.parent()
        powerval = selfval = self.valuation_bigoh(lazylimit=limit_lazyprecision)
        powerprec = selfprec = self._precision
        n = exp
        while n & 1 == 0:
            powerprec = powerval * powerprec + powerprec * powerval
            powerval = powerval * powerval
            n >>= 1
        val = powerval
        precision = powerprec
        n >>= 1
        while n > 0:
            if n & 1 == 1:
                precision = val * powerprec + precision * powerval
                val = val * powerval
            powerprec = powerval * powerprec + powerprec * powerval
            powerval = powerval * powerval
            n >>= 1
        vself = selfval.flat_below()
        if self._group.lazy_mode() is not LazyMode_Disable:
            from precision.misc import ceil_infty
            workprec_shift = (exp-1)*vself
            workprec_transformation = lambda x: max(x - workprec_shift, ceil_infty(x/exp))
            approximation = self._approximation._pow_(exp, workprec_transformation)
        else:
            approximation = self._approximation._pow_(exp)
        return self._new_fast(parent, approximation, precision, cap_precision=True)

    def nrows(self):
        return self.parent().nrows()

    def ncols(self):
        return self.parent().ncols()

    def valuation_bigoh(self,model=None,lazylimit=Infinity):
        from bigoh import FlatBigOh_matrix, JaggedBigOh_matrix, LatticeBigOh_matrix
        valuations = [ self[i,j].valuation_bigoh(lazylimit=lazylimit) for i in range(self.nrows()) for j in range(self.ncols()) ]
        parentprec = self._precision.parent()
        if model is None:
            model = self.precision_model()
        if model is FlatBigOh_matrix:
            valuation = self.parent()._exact_precision
            for prec in valuations:
                valuation += prec
            return FlatBigOh_matrix(parentprec, valuation)
        elif model is JaggedBigOh_matrix or model is ExactBigOh:
            return JaggedBigOh_matrix(parentprec, valuations)
        else:
            raise TypeError("invalid model")

    def _repr_(self):
        if self._approximation.length() >= 0:
            from printing import repr_matrix
            coeffs = [ [ self[i,j] for j in range(self.ncols()) ] for i in range(self.nrows()) ]
            return repr_matrix(coeffs)
        else:
            return Element_inexact._repr_(self)

    def __invert__(self):
        raise NotImplementedError
