#include "/Users/xavier/sage/devel/sage/sage/ext/stdsage.pxi"

from operator import add, mul

from sage.misc.cachefunc import cached_function
from sage.structure.element import Element 
from sage.rings.infinity import Infinity

from bigoh import BigOh
from lazy import LazyApproximation, LazyApproximationWrapper

from grouping.mode import PrecisionMode_Disable, PrecisionMode_Individual, PrecisionMode_Collective
from grouping.mode import LazyMode_Disable, LazyMode_ASAP, LazyMode_Printing, LazyMode_Lazy
from grouping.exception import GroupError

from limits import limit_lazyprecision 
from exception import PrecisionError, ApproximationError


class Element_inexact(Element):
    # for cdef classes, don't forget:
    # cdef object __weakref__

    def __init__(self, parent, x, prec=None, workprec=None, group=None):
        from parent_inexact import Parent_inexact
        if not isinstance(parent, Parent_inexact):
            raise TypeError("parent must be an instance of Parent_inexact: %s" % parent)
        Element.__init__(self,parent)
        base = parent.base_ring()
        baseapprox = base.approximation()
        parentapprox = parent.approximation()

        # determine group
        if group is None:
            from grouping.context import stack
            group = stack[-1]
        else:
            from grouping.group import GroupElements_inexact
            if not isinstance(group, GroupElements_inexact):
                raise TypeError("group not valid")
        self._group = group
        if group.precision_mode() is PrecisionMode_Disable:
            self._precision = None
        else:
            self._precision = parent._exact_precision

        # build approximation
        coeffs = [ ]
        if isinstance(x,list):
            dimension = self.parent().dimension()
            if len(x) > dimension:
                raise ValueError("list too long")
            for i in range(len(x)):
                try:
                    if x[i].parent() is base:
                        coeffs.append(x[i])
                    else:
                        coeffs.append(base(x[i], group=self._group))
                except (TypeError,ValueError,PrecisionError):
                    raise TypeError("unable to create an element from %s" % x)
            approximation = [ c._approximation for c in coeffs ]
            if group.precision_mode() is not PrecisionMode_Disable:
                self._precision = parent.precision()([ c._precision for c in coeffs ])
        elif isinstance(x,dict):
            raise NotImplementedError
        elif isinstance(x,Element_inexact) and x.parent() is parent:
            approximation = x._approximation
            if group.precision_mode() is not PrecisionMode_Disable:
                self._precision = x._precision
        elif isinstance(x,Element) and x.parent() is parentapprox:
            approximation = x
        else:
            raise TypeError("unable to create an element from %s" % x)

        # build precision
        if group.precision_mode() is PrecisionMode_Individual:
            if prec is not None:
                if isinstance(prec,tuple):
                    prec = list(prec)
                else:
                    prec = [ prec ]
                self._precision += parent.precision()(*prec)

        # lazy_mode
        if not isinstance(approximation, (LazyApproximation,LazyApproximationWrapper)):
            if group.lazy_mode() is not LazyMode_Disable:
                approximation = parentapprox._lazy_class(parentapprox,approximation)   # need to split it
            else:
                approximation = parentapprox(approximation)
        if group.lazy_mode() is LazyMode_Disable:
            self._approximation = approximation
        else:
            self._approximation = approximation.new_wrapper()

        # cap precision and set working precision
        if group.precision_mode() is PrecisionMode_Individual:
            self._cap_precision()
            if workprec is None:
                # Bug: for polynomials, to_workprec needs an extra argument
                workprec = self._precision.to_workprec()
        else:
            if workprec is None:
                workprec = Infinity
        self._approximation.update_max_workprec(workprec)
        if group.lazy_mode() is LazyMode_Disable:
            self._approximation = self._approximation.at_workprec()

        # register in group
        group.register_member(self,safemode=True)

        # evaluate
        if group.lazy_mode() is LazyMode_ASAP:
            try:
                self.evaluate()
            except ApproximationError:
                pass

    def _new_fast(self, parent, approximation, precision, cap_precision=True):
        #element = PY_NEW_SAME_TYPE(self)
        cls = self.__class__; element = cls.__new__(cls)
        Element.__init__(element, parent)
        element._group = self._group
        if approximation.is_lazy():
            element._approximation = approximation.new_wrapper()
        else:
            element._approximation = approximation
        element._precision = precision
        if cap_precision:
            element._cap_precision()
        # Bug: for polynomials, to_workprec needs an extra argument
        workprec = precision.to_workprec()
        element._approximation.update_max_workprec(workprec)
        if element._group.lazy_mode() is LazyMode_ASAP:
            try:
                element.evaluate()
            except ApproximationError:
                pass
        return element

    def _cap_precision(self):
        capped_relative = self._group.capped_relative()
        if capped_relative is not Infinity:
            val_bigoh = self.valuation_bigoh(lazylimit=limit_lazyprecision) << capped_relative
            self._precision += val_bigoh

    def __hash__(self):
        # Do NOT use __repr__ here! (since __repr__ evaluates)
        return id(self)
        return hash((self._approximation, self._precision))

    def evaluate(self):
        # Bug: for polynomials, to_workprec needs an extra argument
        self._approximation.at_workprec(self._precision.to_workprec())
        return self

    def group(self):
        return self._group

    def approximation(self, precision=None):
        if precision is None:
            # Bug: for polynomials, to_workprec needs an extra argument
            workprec = self._precision.to_workprec()
        elif isinstance(precision, BigOh) and precision.parent() is self.parent().precision():
            workprec = (self._precision + self.parent().precision()(precision)).to_workprec()
        else:
            workprec = precision
        return self._approximation.at_workprec(workprec)

    def precision(self):
        if self._group.precision_mode() is PrecisionMode_Disable:
             raise PrecisionError("Precision of this element is not tracked")
        return self._precision

    def change_group(self,group=None):
        if group is None:
            from grouping.context import stack
            group = stack[-1]
        else:
            from grouping.group import GroupElements_inexact
            if not isinstance(group, GroupElements_inexact):
                raise TypeError("group not valid")
        return self.__class__(self.parent(), self, group=group)

    def change_precision(self,prec):
        parentprec = self.parent().precision()
        precision = parentprec(prec)
        return self._new_fast(self.parent(), self._approximation, precision, cap_precision=True)

    def add_bigoh(self,prec):
        parentprec = self.parent().precision()
        precision = self._precision + parentprec(prec)
        return self._new_fast(self.parent(), self._approximation, precision, cap_precision=False)

    def precision_model(self):
        if self._precision is None:
            raise PrecisionError("Precision of this element is not tracked")
        else:
            return type(self._precision)

    def working_precision(self):
        if self._group is None and self._workprec is None:
             raise PrecisionError("Working precision of this element is not tracked")
        if self._group is None:
            return self._workprec
        else:
            return self._group.working_precision()

    def is_exact(self):
        if self._group.precision_mode() is PrecisionMode_Disable:
             raise PrecisionError("Precision of this element is not tracked")
        return self._precision.is_exact()

    def is_evaluated(self):
        approximation = self._approximation
        if approximation.is_lazy():
            return approximation.is_evaluated()
        else:
            return True

    def __getitem__(self,i):
        parent = self.parent()
        indices = parent.indices_basis()
        if indices is not None:
            try:
                i = indices.index(i)
            except ValueError:
                raise IndexError("invalid index: %s" % i)
        if i < 0 or i >= parent.dimension():
            raise IndexError("invalid index: %s" % i)
        return self.parent().base_ring()(self._approximation._getitem_by_num(i),prec=self._precision._getitem_by_num(i))

    def __setitem__(self,i,x):
        raise TypeError("Elements are immutable")

    def list(self):
        dimension = self.parent().dimension()
        if dimension is Infinity:
            raise NotImplementedError
        return [ self._getitem_by_num(i) for i in range(dimension) ]

    def valuation_bigoh(self,model=None,lazylimit=Infinity): # here minimal=True
        raise NotImplementedError

    def __str__(self):
        return self._repr_()

    def __repr__(self, truncate=None):
        if self._group.lazy_mode() is LazyMode_Printing:
            try:
                self.evaluate()
            except ApproximationError:
                pass
        s = self._repr_()
        lines = s.split('\n')
        if len(lines) < 2:
            return s
        if truncate is None:
            from window import width
            truncate = width() - 1
        truncated = False
        s = ""
        for line in lines:
            if len(line) > truncate:
                truncated = True
                line = line[:truncate]
            s += line + "\n"
        s = s[:-1]
        if truncated:
            from sage.misc.sageinspect import sage_getvariablename
            first_line = "Output truncated"
            name = sage_getvariablename(self)
            if isinstance(name, str):
                first_line += " (type 'print %s.__repr__(truncate=Infinity)' to see the full output)" % name
            return first_line + ":\n" + s
        else:
            return s

    def _repr_(self):
        from printing import repr_operator
        approximation = "%s" % self._approximation
        if self.is_exact():
            return approximation
        else:
            precision = "%s" % self._precision
            return repr_operator(" + ", [approximation, precision])

    def is_zero(self,sure=False,lazylimit=Infinity):
        if sure:
            return self.is_exact() and self._approximation.is_zero()
        else:
            bigoh = self._precision
            #if self.is_lazy():
            bigoh += self.parent().precision()(lazylimit)
            return self.valuation_bigoh(lazylimit=lazylimit) <= bigoh
 
    #def __richcmp__(self):


class ModuleElement_inexact(Element_inexact):
    def __add__(self,other):
        parent = self.parent()
        if parent == other.parent():
            group = self.group()
            if other.group() != group:
                raise GroupError("Can add only two elements in the same group")
            return self._add_(other)
        else:
            from sage.structure.element import get_coercion_model
            return get_coercion_model().bin_op(self,other,add)

    def _add_(self,other):
        parent = self.parent()
        precision = self._precision + other._precision
        if self._group.lazy_mode() is not LazyMode_Disable:
            expected_valuation = self.valuation_bigoh(lazylimit=limit_lazyprecision) + other.valuation_bigoh(lazylimit=limit_lazyprecision)
            approximation = self._approximation._add_(other._approximation, expected_valuation = expected_valuation.flat_below())
        else:
            approximation = self._approximation._add_(other._approximation)
        return self._new_fast(parent, approximation, precision)

    def __neg__(self):
        parent = self.parent()
        if self._group.lazy_mode() is not LazyMode_Disable:
            approximation = self._approximation.__neg__()
        else:
            approximation = self._approximation.__neg__()
        return self._new_fast(parent, approximation, self._precision, cap_precision=False)

    def __sub__(self, other):
        return self + (-other)

    def _acted_upon_(self, x, self_on_left):
        coerce_from_base = False
        base = self.base_ring()
        if x.parent() is base:
            coerce_from_base = True
        else:
            try:
                x = base(x)
                coerce_from_base = True
            except:
                pass
        if coerce_from_base:
            if x.group() != self.group():
                raise GroupError("Can multiply only two elements in the same group")
            if self_on_left:
                return self._lmul_(x)
            else:
                return self._rmul_(x)

    def _rmul_(self,coeff):  # coeff * self
        parent = self.parent()
        valself = self.valuation_bigoh(lazylimit=limit_lazyprecision)
        valcoeff = coeff.valuation_bigoh(lazylimit=limit_lazyprecision)
        precision = self._precision._rmul_(valcoeff) + valself._rmul_(coeff._precision)
        if self._group.lazy_mode() is not LazyMode_Disable:
            expected_valuation = valself._rmul_(valcoeff)
            approximation = self._approximation._rmul_(coeff._approximation, expected_valuation = expected_valuation.flat_below(), workprec_shifts = [ -valcoeff.flat_below(), -valself.flat_below() ])
        else:
            approximation = self._approximation._rmul_(coeff._approximation)
        return self._new_fast(parent, approximation, precision, cap_precision=False)

    def _lmul_(self,coeff):  # self * coeff
        parent = self.parent()
        valself = self.valuation_bigoh(lazylimit=limit_lazyprecision)
        valcoeff = coeff.valuation_bigoh(lazylimit=limit_lazyprecision)
        precision = self._precision._lmul_(valcoeff) + valself._lmul_(coeff._precision)
        if self._group.lazy_mode() is not LazyMode_Disable:
            expected_valuation = valself._lmul_(valcoeff)
            approximation = self._approximation._lmul_(coeff._approximation, expected_valuation = expected_valuation.flat_below(), workprec_shifts = [ -valcoeff.flat_below(), -valself.flat_below() ])
        else:
            approximation = self._approximation._lmul_(coeff._approximation)
        return self._new_fast(parent, approximation, precision, cap_precision=False)


class RingElement_inexact(ModuleElement_inexact):
    def __mul__(self,other):
        parent = self.parent()
        if parent == other.parent():
            group = self.group()
            if other.group() != group:
                raise GroupError("Can multiply only two elements in the same group")
            return self._mul_(other)
        else:
            from sage.structure.element import get_coercion_model
            return get_coercion_model().bin_op(self,other,mul)

    def _mul_(self,other):
        parent = self.parent()
        selfval = self.valuation_bigoh(lazylimit=limit_lazyprecision)
        otherval = other.valuation_bigoh(lazylimit=limit_lazyprecision)
        precision = self._precision * otherval + selfval * other._precision
        if self._group.lazy_mode() is not LazyMode_Disable:
            expected_valuation = selfval * otherval
            approximation = self._approximation._mul_(other._approximation, expected_valuation = expected_valuation.flat_below(), workprec_shifts = [ -otherval.flat_below(), -selfval.flat_below() ])
        else:
            approximation = self._approximation._mul_(other._approximation)
        return self._new_fast(parent, approximation, precision)

    def __pow__(self,exp,ignored=None):
        if exp < 0:
            return (~self).__pow__(-exp)
        elif exp == 0:
            return self.parent().one()
        elif exp == 1:
            return self
        parent = self.parent()
        if parent.is_commutative():
            selfval = self.valuation_bigoh(lazylimit=limit_lazyprecision)
            precision = selfval**(exp-1) * self._precision
            val = self.base_ring()(exp).valuation()
            precision <<= val
            precision += self._precision ** exp
            vself = selfval.flat_below()
        else:
            powerval = selfval
            powerprec = selfprec = self._precision
            n = exp
            while n & 2 == 0:
                powerprec = powerval * powerprec + powerprec * powerval
                powerval = powerval * powerval
                n >>= 1
            val = powerval
            precision = powerprec
            n >>= 1
            while n > 0:
                if n & 2 == 1:
                    precision = val * powerprec + precision * powerval
                    val = val * powerval
                powerprec = powerval * powerprec + powerprec * powerval
                powerval = powerval * powerval
                n >>= 1
            vself = selfval.flat_below()
            val = 0
        if self._group.lazy_mode() is not LazyMode_Disable:
            from sage.functions.other import ceil
            workprec_shift = (exp-1)*vself + val
            def f(x):
                if x == Infinity: return Infinity
                else: return max(x - workprec_shift, ceil(x/exp))
            approximation = self._approximation._pow_(exp, workprec_transformation=f, expected_valuation=vself*exp)
        else:
            approximation = self._approximation._pow_(exp)
        return self._new_fast(parent, approximation, precision, cap_precision=True)

    def __invert__(self):
        raise NotImplementedError

    def __div__(self,other):
        return self * (~other)

    def is_one(self,sure=False):
        return (self - self.parent().one()).is_zero(sure=sure)
