from sage.rings.infinity import Infinity

from precision.lazy import lazy_method
from precision.lazy import LazyApproximation
from approximation import MatrixApprox


class LazyApproximation_matrix(LazyApproximation):
    def __init__(self, parent, function, starting_workprec=None, **kwargs):
        if isinstance(function, list):
            nrows = parent.nrows()
            ncols = parent.ncols()
            if len(function) > nrows*ncols:
                raise ValueError("list too long")
            coeffs = function + ([0] * (nrows*ncols-len(function)))
            coeffs = [ coeffs[i*ncols:(i+1)*ncols] for i in range(nrows) ]
            from printing import repr_matrix
            self.repr = lambda: repr_matrix(coeffs)
        LazyApproximation.__init__(self, parent, function, starting_workprec, **kwargs)

    #def valuation_entries(self,lazylimit=Infinity):
    #    workprec = self._current_workprec
    #    shift = 1
    #    nrows = self.parent().nrows()
    #    ncols = self.parent().ncols()
    #    while True:
    #        approximation = self.at_workprec(workprec)
    #        vals = [ ]
    #        not_enough_prec = False
    #        for i in range(nrows):
    #            for j in range(ncols):
    #                val = approximation[i,j].valuation()
    #                if val < workprec:
    #                    vals.append(val)
    #                else:
    #                    not_enough_prec = True
    #        if not not_enough_prec:
    #            return vals
    #        if workprec >= lazylimit:
    #            return [ min(workprec, val) for val in vals ]
    #        workprec += shift
    #        shift *= 2
