from operator import add, mul

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element import Element
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.infinity import Infinity

from parent_precision import ParentBigOh
from coerce import coerce_bigoh
from exception import PrecisionError


class BigOh(Element):
    def __init__(self,parent,*args,**kwargs):
        if not isinstance(parent,ParentBigOh):
            raise TypeError("parent must be an instance of ParentBigOh")
        Element.__init__(self,parent)
        self._baseprec = parent.base_precision()
        if self._baseprec is None:
            self._base_exactprec = None
        else:
            if self._baseprec is parent and isinstance(self,ExactBigOh):
                self._base_exactprec = self
            else:
                self._base_exactprec = ExactBigOh(self._baseprec)
        if kwargs.has_key('convert'):
            convert = kwargs['convert']
            del kwargs['convert']
        else:
            convert = (len(args) > 0 and isinstance(args[0],BigOh) and args[0].parent() is parent)
        if convert:
            self._convert_(parent,*args,**kwargs)
        else:
            self._init_(parent,*args,**kwargs)

    def _init_(self,parent,*args,**kwargs):
        pass

    def _convert_(self,parent,prec):
        raise NotImplementedError("conversion %s -> %s not implemented" % (prec._repr_model(),self._repr_model()))

    def parenthesis_level(self):
        return 3

    def _repr_model(self):
        s = str(self.__class__)
        return s[s.rfind(".")+1:-2]

    def __add__(self,other):
        if self.parent() != other.parent(): 
            # should implement coercion here
            raise NotImplementedError
        if self.is_exact():
            return other
        if other.is_exact():
            return self
        if self.__class__ is other.__class__:
            return self._add_(other)
        return coerce_bigoh(self,other,add)

    def _add_(self,other):
        raise NotImplementedError

    def intersection(self,other):
        if self.parent() != other.parent(): 
            # should implement coercion here
            raise NotImplementedError
        if self.is_exact():
            return self
        if other.is_exact():
            return other
        if self.__class__ is other.__class__:
            return self._intersection(other)
        return coerce_bigoh(self,other,self.__class__.intersection)

    def _intersection(self,other):
        raise NotImplementedError

    def __mul__(self,other):
        if self.parent() != other.parent():
            # should implement coercion here
            raise NotImplementedError
        if self.is_exact():
            return self
        if other.is_exact():
            return other
        if self.__class__ is other.__class__:
            return self._mul_(other)
        return coerce_bigoh(self,other,mul)

    def _mul_(self,other):
        raise NotImplementedError

    def __pow__(self,exp,ignored=None):
        if exp > 0:
            power = self
            while exp & 1 == 0:
                power *= power
                exp >>= 1
            exp >>= 1
            res = power
            while exp > 0:
                power *= power
                if exp & 1 == 1:
                    res *= power
                exp >>= 1
            return res
        else:
            raise ValueError("The exponent must be a positive integer")

    def __lshift__(self,i):
        raise NotImplementedError
    def __rshift__(self,i):
        raise NotImplementedError

    def is_exact(self):
        return self.flat_below() is Infinity

    def first(self):
        return 0

    def last(self):
        return self.parent().dimension()

    def flat_below(self):
        raise NotImplementedError

    def to_workprec(self):
        raise NotImplementedError

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
        return self._getitem_by_num(i)

    def _getitem_by_num(self,i):
        raise NotImplementedError

    def __repr__(self,model=False):
        if model:
            return self._repr_model() + ": " + self._repr_(model=True)
        else:
            return self._repr_(model=False)

    def _repr_(self,model):
        return "Generic BigOh in %s" % self._ambiant_space

    def contains(self,other):
        raise NotImplementedError

    def equals(self,other):
        return self.__repr__() == other.__repr__()

    def __eq__(self,other):
        if self.parent() != other.parent():
            raise NotImplementedError("Parents differ in comparison of bigoh's")
        if self.is_exact():
            return other.is_exact()
        if other.is_exact():
            return self.is_exact()
        if self.__class__ is not other.__class__:
            raise NotImplementedError("Models differ in comparison of bigoh's")
        return self.equals(other)

    def __ne__(self,other):
        return not (self == other)

    def __ge__(self,other):
        if self.parent() != other.parent():
            raise NotImplementedError("Parents differ in comparison of bigoh's")
        if self.is_exact():
            return other.is_exact()
        if other.is_exact():
            return True
        if self.__class__ is not other.__class__:
            raise NotImplementedError("Models differ in comparison of bigoh's")
        return self.contains(other)

    def __gt__(self,other):
        return self >= other and not self.equals(other)

    def __le__(self,other):
        return other >= self

    def __lt__(self,other):
        return other >= self and not self.equals(other)


class ExactBigOh(BigOh,UniqueRepresentation):
    def flat_below(self):
        return Infinity

    def to_workprec(self):
        return Infinity

    def _add_(self,other):
        return other

    def _intersection(self,other):
        return self

    def _mul_(self,other):
        return self

    def _lmul_(self,coeff):
        return self

    def _rmul_(self,coeff):
        return self

    def __lshift__(self,i):
        return self
    def __rshift__(self,i):
        return self

    def sqrt(self,exp=2):
        return self

    def is_exact(self):
        return True

    def _repr_(self,model):
        return "Exact"

    def first(self):
        return self.parent().dimension()

    def last(self):
        return 0

    def _getitem_by_num(self,i):
        return ExactBigOh(ParentBigOh(self.parent().base_ring()))


class ValuationBigOh(BigOh):
    def _init_(self,parent,prec):
        if prec is Infinity:
            self._precision = Infinity
        else:
            self._precision = ZZ(prec)

    def _convert_(self,parent,prec):
        self._precision = prec.flat_below()

    def flat_below(self):
        return self._precision

    def to_workprec(self):
        return self._precision

    def _add_(self,other):
        return ValuationBigOh(self.parent(), min(self._precision,other._precision))

    def _intersection(self,other):
        return ValuationBigOh(self.parent(), max(self._precision,other._precision))

    def _mul_(self,other):
        return ValuationBigOh(self.parent(), self._precision + other._precision)

    def _lmul_(self,coeff):
        return self._mul_(coeff)

    def _rmul_(self,coeff):
        return self._mul_(coeff)

    def __lshift__(self,i):
        return ValuationBigOh(self.parent(), self._precision + i)
    def __rshift__(self,i):
        return ValuationBigOh(self.parent(), self._precision - i)

    def sqrt(self,exp=2):
        from sage.functions.other import ceil
        precision = ceil(self._precision / exp)
        return ValuationBigOh(self.parent(), precision)

    def is_exact(self):
        return self._precision is Infinity

    def _getitem_by_num(self,i):
        return self

    def equals(self,other):
        return self._precision == other._precision

    def contains(self,other):
        return self._precision <= other._precision

    def _repr_(self,model):
        if self._precision == Infinity:
            return "Exact"
        elif self._precision == 0:
            return "O(1)"
        elif self._precision == 1:
            return "O(%s)" % self.parent().uniformizer_name()
        else:
            return "O(%s^%s)" % (self.parent().uniformizer_name(), self._precision)


class FlatBigOh(BigOh):
    def _init_(self,parent,prec):
        self._precision = parent.base_precision()(prec)

    def _convert_(self,parent,prec):
        self._precision = ExactBigOh(parent.base_precision())
        for i in range(self.parent.dimension()):
            self._precision += prec._getitem_by_num(i)

    def flat_below(self):
        return self._precision.flat_below()

    def to_workprec(self):
        return self._precision.to_workprec()

    def _add_(self,other):
        return self.__class__(self.parent(), self._precision + other._precision)

    def _intersection(self,other):
        return self.__class__(self.parent(), self._precision.intersection(other._precision))

    def _mul_(self,other):
        return self.__class__(self.parent(), self._precision * other._precision)

    def _lmul_(self,coeff):
        return self.__class__(self.parent(), self._precision * coeff)

    def _rmul_(self,coeff):
        return self.__class__(self.parent(), coeff * self._precision)

    def __lshift__(self,i):
        return self.__class__(self.parent(), self._precision.__lshift__(i))
    def __rshift__(self,i):
        return self.__class__(self.parent(), self._precision.__rshift__(i))

    def is_exact(self):
        return self._precision.is_exact()

    def _getitem_by_num(self,i):
        return self._precision

    def equals(self,other):
        return self._precision == other._precision

    def contains(self,other):
        return self._precision >= other._precision

    def _repr_(self,model):
        return self._precision._repr_(model=model)


class JaggedBigOh(BigOh):
    def _init_(self,parent,prec):
        dimension = self.parent().dimension()
        if isinstance(prec,list):
            if len(prec) > dimension:
                raise ValueError("list too long")
            self._precvector = [ parent.base_precision()(p) for p in prec ]
        elif isinstance(prec,dict):
            indices = self.parent().indices_basis()
            max = 0
            for key in prec:
                if indices is not None:
                    try:
                        i = indices.index(key)
                    except IndexError:
                        raise IndexError("invalid index: %s" % key)
                else:
                    i = key
                    if i < 0 or i >= dimension:
                        raise IndexError("invalid index: %s" % i)
                if i > max:
                    self._precvector += (max-i-1) * [self._base_exactprec]
                    self._precvector.append(prec[key])
                    max = i
                else:
                    self._precvector[i] = prec[key]
        elif dimension < Infinity:
            base = ParentBigOh(parent.base_ring())
            self._precvector = dimension * [base(prec)]
        else:
            raise TypeError("prec must be a list or a dictionary")
        self._normalize()

    def _convert_(self,parent,prec):
        self._precvector = [ prec._getitem_by_num(i) for i in range(prec.last()) ]

    def _normalize(self):
        for i in range(len(self._precvector)-1,-1,-1):
            if self._precvector[i].is_exact():
                del self._precvector[i]
            else:
                break

    def _add_(self,other):
        l1 = len(self._precvector)
        l2 = len(other._precvector)
        if l1 <= l2:
            prec = [ self._precvector[i] + other._precvector[i] for i in range(l1) ] + other._precvector[l1:]
        else:
            prec = [ self._precvector[i] + other._precvector[i] for i in range(l2) ] + self._precvector[l2:]
        return self.__class__(self.parent(), prec)

    def _lmul_(self,coeff):
        return self.__class__(self.parent(), [ prec * coeff for prec in self._precvector ])

    def _rmul_(self,coeff):
        return self.__class__(self.parent(), [ coeff * prec for prec in self._precvector ])

    def _intersection(self,other):
        l = min(len(self._precvector), len(other._precvector))
        prec = [ self._precvector[i].intersection(other._precvector[i]) for i in range(l) ]
        return self.__class__(self.parent(), prec)

    def to_workprec(self,last=None):
        dimension = self.parent().dimension()
        if last is None:
            if dimension is Infinity:
                last = self.last()
            else:
                last = dimension
        length = len(self._precvector)
        if last < length:
            raise ValueError("still inexact BigOh's after last")  # is it really an error?
        if length < last:
            return Infinity
        if length > 0:
            maximum = -Infinity
            for p in self._precvector:
                if p.is_exact():
                    maximum = Infinity
                    break
                v = p.to_workprec()
                if v > maximum: maximum = v
            return maximum
        else:
            return Infinity

    def flat_below(self):
        if len(self._precvector) > 0:
            return min([ p.flat_below() for p in self._precvector ])
        else:
            return -Infinity

    def __lshift__(self,i):
        prec = [ p.__lshift__(i) for p in self._precvector ]
        return self.__class__(self.parent(), prec)
    def __rshift__(self,i):
        prec = [ p.__rshift__(i) for p in self._precvector ]
        return self.__class__(self.parent(), prec)

    def is_exact(self):
        return len(self._precvector) == 0

    def _getitem_by_num(self,i):
        if i < len(self._precvector):
            return self._precvector[i]
        else:
            return self._base_exactprec

    def first(self):
        l = len(self._precvector)
        for i in range(l):
            if not self._precvector[i].is_exact():
                return i
        return l

    def last(self):
        return len(self._precvector)

    def equals(self,other):
        length = len(self._precvector)
        if len(other._precvector) != length:
            return False
        for i in range(length):
            if self._precvector[i] != other._precvector[i]:
                return False
        return True

    def contains(self,other):
        length = len(self._precvector)
        if len(other._precvector) > length:
            return False
        for i in range(length):
            if not (self._precvector >= other._precvector[i]):
                return False
        return True

    def _repr_(self,model=False):
        if len(self._precvector) == 0:
            return "Exact"
        s = ""
        if model:
            indent = "    "
        else:
            indent = ""
        indices = self.parent().indices_basis()
        if indices is None:
            for p in self._precvector:
                s += "%s, " % p
            if self.parent().dimension() < Infinity:
                s += "Exact, " * (self._indices - len(self._precvector))
                s = s[:-2]
            else:
                s += "Exact, ..."
        else:
            for i in range(len(self._precvector)):
                key = indices[i]
                prec = self._precvector[i]
                if not prec.is_exact():
                    s += "\n%s%s: %s" % (indent,key,self._precvector[i])
            if not model: s = s[1:]
        return s


class FlatIntervalBigOh(BigOh):
    def _init_(self,parent,prec,last=None,first=0):
        dimension = self.parent().dimension()
        if first < 0:
            raise ValueError("first must be a nonnegative integer")
        if last is None:
            last = dimension
        elif last > dimension:
            raise ValueError("last must be less than or equal to the dimension (=%s)" % dimension)
        self._precision = parent.base_precision()(prec)
        if self._precision.is_exact():
            self._first = dimension
            self._last = 0
        else:
            self._first = first
            self._last = last

    def _convert_(self,parent,prec):
        self._first = prec.first()
        self._last = prec.last()
        if isinstance(prec,FlatBigOh):
            self._precision = prec[0]
        else:
            self._precision = self._base_exactprec
            for i in range(self._first, self._last):
                self._precision += prec[i]

    def flat_below(self):
        if self._last > self._first:
            return self._precision.flat_below()
        else:
            return Infinity

    def to_workprec(self,last=None):
        dimension = self.parent().dimension()
        if last is None:
            if dimension is Infinity:
                last = self._last
            else:
                last = dimension
        if last < self._last:
            raise ValueError("still inexact BigOh's after last")  # is it really an error?
        if last <= self._first:
            return Infinity
        elif self._first == 0:
            return self._precision.to_workprec()
        else:
            return Infinity

    def _add_(self,other):
        return FlatIntervalBigOh(self.parent(), self._precision + other._precision, max(self._last, other._last), min(self._first, other._first))

    def _lmul_(self,coeff):
        return FlatIntervalBigOh(self.parent(), self._precision * coeff, self._last, self._first)

    def _rmul_(self,coeff):
        return FlatIntervalBigOh(self.parent(), coeff * self._precision, self._last, self._first)

    def _intersection(self,other):
        return FlatIntervalBigOh(self.parent(), self._precision.intersection(other._precision), min(self._last, other._last), max(self._first, other._first))

    def __lshift__(self,i):
        return FlatIntervalBigOh(self.parent(), self._precision.__lshift__(i), self._last, self._first)
    def __rshift__(self,i):
        return FlatIntervalBigOh(self.parent(), self._precision.__rshift__(i), self._last, self._first)

    def is_exact(self):
        if self._first >= self._last:
            return True
        return self._precision.is_exact()

    def _getitem_by_num(self,i):
        if i >= self._first and i <= self._last:
            return self._precision
        else:
            return self._base_exactprec

    def first(self):
        return self._first

    def last(self):
        return self._last

    def equals(self,other):
        # Assume self and other are not exact (checked in __eq__)
        return self._first == other._first and self._last == other._last and self._precision == other._precision

    def contains(self,other):
        # Assume self and other are not exact (checked in __ge__)
        return self._first <= other._first and self._last >= other._last and self._precision >= other._precision

    def _repr_(self,model):
        if self._precision is self._base_exactprec or self._last <= self._first:
            return "Exact"
        else:
            return "%s on interval [%s,%s]" % (self._precision, self._first, self._last-1)


class LatticeBigOh(BigOh):
    pass
