from sage.rings.rational_field import QQ
from sage.rings.infinity import Infinity

from sage.structure.element import RingElement

from precision.element_inexact import RingElement_inexact
from approximation import QQpApprox_element
from precision.bigoh import ExactBigOh, ValuationBigOh
from lazy import LazyApproximation_padics
from precision.grouping.mode import LazyMode_Disable

from precision.limits import limit_lazyprecision
from precision.exception import PrecisionError, ApproximationError


class QQp_element(RingElement_inexact):
    def __init__(self, parent, x, prec=None, group=None, *args, **kwargs):
        from ZZp_element import ZZp_element
        approxring = parent.approximation()
        if isinstance(x, ZZp_element):
            approximation = x._approximation
            prec = x._precision.flat_below()
        elif x in QQ:
            approximation = QQpApprox_element(approxring, x, 0)
        else:
            approximation = x
        RingElement_inexact.__init__(self, parent, approximation, prec, group=group, *args, **kwargs)

    def parenthesis_level(self):
        if self.is_exact():
            try:
                return self._approximation.parenthesis_level()
            except (AttributeError, TypeError):
                return 3
        else:
            return 0

    def lift_to_precision(self,prec):
        return QQp_element(self.parent(), self._approximation, prec, group=self._group)

    def _repr_(self):
        prec = self._precision.flat_below()
        if prec is Infinity:
            return "%s" % self._approximation
        else:
            c = self._approximation
            if c.is_evaluated_zero():
                s = ""
            else:
                if c.is_evaluated():
                    c = c.at_workprec()
                s = "%s + " % c
            if prec == 0:
                return s + "O(1)"
            elif prec == 1:
                return s + "O(%s)" % (self.parent().uniformizer_name())
            else:
                return s + "O(%s^%s)" % (self.parent().uniformizer_name(), prec)

    def valuation(self,minimal=False,lazylimit=Infinity):
        v = self._approximation.valuation(lazylimit=lazylimit)
        prec = self._precision.flat_below()
        if minimal:
            v = min(v,prec)
        else:
            if v >= prec: v = Infinity
        return v

    def valuation_bigoh(self,model=None,lazylimit=Infinity):
        prec = self._precision.flat_below()
        v = min(prec, self._approximation.valuation(lazylimit=lazylimit))
        if model is ValuationBigOh or model is ExactBigOh or model is None:
            return self._precision.parent()(v)
        else:
            raise TypeError("invalid model")

    def is_unit(self):
        return self.valuation() == 0

    def __invert__(self):
        val = self.valuation()
        if val is Infinity:
            raise ZeroDivisionError
        precision = self._precision >> (2*val)
        parent = self.parent()
        if self._group.lazy_mode() is not LazyMode_Disable:
            approximation = self._approximation._invert_(workprec_shift=2*val, expected_valuation=-val)
        else:
            approximation = self._approximation._invert_()
        return self._new_fast(parent, approximation, precision, cap_precision=False)

    def log(self,lazylimit=None):
        parent = self.parent()
        x = self - parent.one()
        if lazylimit is None:
            lazylimit = max(limit_lazyprecision, 1)
        val = x.valuation(lazylimit=lazylimit,minimal=True)
        if val is Infinity:
            return self.parent()(0)
        p = parent._p
        if val == 0:
            raise NotImplementedError("The argument must be congruent to 1 modulo %s" % p)
        precision = self._precision
        if self._group.lazy_mode() is not LazyMode_Disable:
            from sage.functions.log import log
            from sage.functions.other import floor
            lazy_options = {
                'expected_valuation': val
            }
            approximation = self._approximation.log(lazy_options=lazy_options)
        else:
            approximation = self._approximation.log(workprec=precision.to_workprec())
        return self._new_fast(parent, approximation, precision, cap_precision=False)

    def exp(self,lazylimit=None):
        parent = self.parent()
        if lazylimit is None:
            lazylimit = max(limit_lazyprecision, 2)
        val = self.valuation(lazylimit=lazylimit,minimal=True)
        if val is Infinity:
            return self.parent()(1)
        p = parent._p
        if p > 2 and val == 0:
            raise ValueError("The argument must be congruent to 0 modulo %s" % p)
        if p == 2 and val <= 1:
            raise ValueError("The argument must be congruent to 0 modulo 4")
        precision = self._precision
        if self._group.lazy_mode() is not LazyMode_Disable:
            from sage.functions.other import ceil
            lazy_options = {
                'expected_valuation': 0
            }
            approximation = self._approximation.exp(lazy_options=lazy_options)
        else:
            approximation = self._approximation.exp(workprec=precision.to_workprec())
        return self._new_fast(parent, approximation, precision, cap_precision=False)

    def teichmuller(self,prec=None):
        parent = self.parent()
        val = self.valuation(lazylimit=1)
        if val < 0:
            raise ValueError("The argument must lie in Z_%s" % parent._p)
        if self._precision.flat_below() <= 0:
            raise PrecisionError("not enough precision")
        if prec is None:
            prec = self._group.capped_relative()
        precision = parent.precision()(prec)
        residue = self.approximation(1)
        if residue.is_zero():
            return self.parent().zero()
        if residue.is_one():
            return self.parent().one()
        if self._group.lazy_mode() is not LazyMode_Disable:
            lazy_options = {
                'expected_valuation': 0,
                'workprec_transformations': [ lambda x: 1 ],
                'repr': lambda: "[%s]" % residue
            }
            approximation = self._approximation.teichmuller(lazy_options=lazy_options)
        else:
            approximation = self._approximation.teichmuller(workprec=precision.to_workprec())
        return self._new_fast(parent, approximation, precision, cap_precision=True)
