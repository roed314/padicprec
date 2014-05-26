from sage.rings.infinity import Infinity
from precision.lazy import LazyApproximation
from precision.lazy import lazy_method
from approximation import QQpApprox_element


class LazyApproximation_padics(LazyApproximation):
    def __init__(self, *args, **kwargs):
        LazyApproximation.__init__(self, *args, **kwargs)
        self._length = -1
        if not hasattr(self, 'expected_valuation'):
            self.expected_valuation = 0

    def valuation(self,lazylimit=Infinity):
        if self._current_workprec is not None:
            workprec = self._current_workprec + 1
        else:
            workprec = 0
        if hasattr(self, 'expected_valuation'):
            workprec = max(workprec, self.expected_valuation + 1)
        shift = 1
        while True:
            approximation = self.at_workprec(workprec)
            val = approximation.valuation()
            if val < workprec:
                return val
            if workprec >= lazylimit and val >= lazylimit:
                return min(workprec, val)
            workprec += shift
            shift *= 2

    def _invert_(self,workprec_shift,expected_valuation=0):
        parent = self.parent()
        def lazy_invert(workprec):
            return self.at_workprec(workprec+workprec_shift).__invert__()
        return self.__class__(parent, lazy_invert, args=[self], expected_valuation=expected_valuation)

    log = lazy_method(QQpApprox_element.log)
    exp = lazy_method(QQpApprox_element.exp)
    teichmuller = lazy_method(QQpApprox_element.teichmuller)
