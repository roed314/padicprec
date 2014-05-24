from sage.rings.rational_field import QQ
from sage.rings.infinity import Infinity

from sage.structure.element import RingElement

from precision.element_inexact import RingElement_inexact
from approximation import QQpApprox_element
from precision.bigoh import ExactBigOh, ValuationBigOh
from QQp_element import QQp_element
from lazy import LazyApproximation_padics

from precision.exception import ApproximationError


class ZZp_element(QQp_element):
    def __init__(self, parent, x, prec=None, *args, **kwargs):
        approxring = parent.approximation()
        if isinstance(x, QQpApprox_element) or isinstance(x, LazyApproximation_padics):
            approximation = x
        elif x in QQ:
            if x.valuation(parent._p) < 0:
                raise ValueError("%s is not in Z_%s" % (x, parent._p))
            approximation = QQpApprox_element(approxring, x, 0)
        else:
            raise TypeError
        RingElement_inexact.__init__(self, parent, approximation, prec, *args, **kwargs)

    def is_unit(self):
        return self.valuation() == 0

    def __invert__(self):
        if self._group is not None:
            raise NotImplementedError
        parent = self.parent()
        if self.valuation() != 0:
            raise ValueError("%s is not invertible in Z_%s" % (self, parent._p))
        precision = self._precision
        def lazy_invert(prec):
            return ~self.approximation(prec)
        default_precision = ValuationBigOh(parent.precision(), 1)
        approximation = parent._lazy_class(parent.approximation(), lazy_invert, default_precision, args=[self])
        return self._new_fast(self.parent(), approximation, precision)
