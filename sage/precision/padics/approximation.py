from sage.structure.unique_representation import UniqueRepresentation

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.infinity import Infinity

from sage.rings.ring import CommutativeRing
from sage.structure.element import RingElement
from precision.approximation import Approximation

from precision.exception import ApproximationError


class QQpApprox(CommutativeRing,UniqueRepresentation):
    def __init__(self, p):
        if not p.is_prime():
            raise ValueError("p must be a prime number")
        self._p = p
        CommutativeRing.__init__(self, self)
        from lazy import LazyApproximation_padics
        self._lazy_class = LazyApproximation_padics
        self._zero = self(0)
        self._lazy_zero = self._lazy_class(self, self._zero)
        self._lazy_zero._length = 0

    def _element_constructor_(self,x,*args,**kwargs):
        if isinstance(x,int): x = ZZ(x)
        return QQpApprox_element(self, x, *args, **kwargs)

    def __repr__(self):
        return "Approximation of %s-adic Ring/Field" % self._p

    def uniformizer(self,exp=1):
        return self._p ** exp

    def dimension(self):
        return 1


class QQpApprox_element(RingElement,Approximation):  # maybe should create a ApproximatedRingElement
    def __init__(self,parent,x,val=0,normalized=False):
        RingElement.__init__(self,parent)
        Approximation.__init__(self,parent)
        self._x = QQ(x)
        self._val = val
        self._normalized = normalized

    def _normalize(self):
        if not self._normalized:
            if self._x == 0:
                self._val = 0
            else:
                powers = self.parent().uniformizer
                v = self._x.valuation(powers())
                self._val += v
                if v > 0:
                    self._x /= powers(v)
                elif v < 0:
                    self._x *= powers(-v)
            self._normalized = True

    def parenthesis_level(self):
        if self._x.denominator() == 1:
            return 3
        else:
            return 1

    def _getitem_by_num(self,i):
        return self

    def is_zero(self):
        return self._x == 0

    def _add_(self,other,**kwargs):
        parent = self.parent()
        selfval = self._val
        otherval = other._val
        if selfval == otherval:
            return QQpApprox_element(parent, self._x + other._x, selfval)
        if selfval > otherval:
            return QQpApprox_element(parent, self._x * parent.uniformizer(selfval-otherval) + other._x, otherval)
        else:
            return QQpApprox_element(self.parent(), self._x + other._x * parent.uniformizer(otherval-selfval), selfval)

    def __neg__(self,**kwargs):
        return QQpApprox_element(self.parent(), -self._x, self._val)

    def _sub_(self,other,**kwargs):
        parent = self.parent()
        selfval = self._val
        otherval = other._val
        if selfval == otherval:
            return QQpApprox_element(parent, self._x - other._x, selfval)
        if selfval > otherval:
            return QQpApprox_element(parent, self._x * parent.uniformizer(selfval-otherval) - other._x, otherval)
        else:
            return QQpApprox_element(self.parent(), self._x - other._x * parent.uniformizer(otherval-selfval), selfval)

    def _mul_(self,other,**kwargs):
        return QQpApprox_element(self.parent(), self._x * other._x, self._val + other._val, self._normalized and other._normalized)

    def _rmul_(self,other,**kwargs):
        return self._mul_(other,**kwargs)

    def _lmul_(self,other,**kwargs):
        return self._mul_(other,**kwargs)

    def __invert__(self,**kwargs):
        return QQpApprox_element(self.parent(), ~self._x, -self._val, self._normalized)

    def _div_(self,other,**kwargs):
        return QQpApprox_element(self.parent(), self._x / other._x, self._val - other._val, self._normalized and other._normalized)

    def valuation(self,lazylimit=Infinity):
        self._normalize()
        if self._x == 0:
            return Infinity
        else:
            return self._val

    def _repr_(self):
        self._normalize()
        if self._x == 0:
            return "0"
        elif self._val == 0:
            return "%s" % self._x
        else:
            return "%s*%s^%s" % (self._x, self.parent()._p, self._val)

    def __pow__(self,exp,**kwargs):
        return QQpApprox_element(self.parent(), self._x ** exp, self._val * exp, self._normalized)

    def __cmp__(self,other):
        v = self._val - other._val
        return cmp(self._x, other._x * self.parent().uniformizer(v))

    def truncate(self,workprec):
        if workprec == Infinity:
            return self
        parent = self.parent()
        self._normalize()
        pow = workprec - self._val
        if pow <= 0:
            return QQpApprox_element(parent, 0, 0, True)
        modulo = parent.uniformizer(pow)
        num = self._x.numerator() % modulo
        denom = self._x.denominator() % modulo
        if denom != 1:
            _,inv,_ = denom.xgcd(modulo)
            num = (num * inv) % modulo
        return QQpApprox_element(parent, num, self._val, True)

    def log(self,workprec=Infinity):
        from sage.functions.log import log
        from sage.functions.other import floor
        if workprec is Infinity:
            raise ApproximationError("unable to compute log to infinite precision")
        parent = self.parent()
        pow = parent(-1)
        res = parent(0)
        t = parent(1) - self
        iter = workprec + floor(log(workprec)/log(parent._p)) + 1
        for i in range(1,iter):
            pow *= t
            res += pow / parent(i)
            res = res.truncate(workprec)
        return res

    def exp(self,workprec=Infinity):
        from sage.functions.other import ceil
        if workprec is Infinity:
            raise ApproximationError("unable to compute exp to infinite precision")
        parent = self.parent()
        pow = parent(1)
        res = parent(1)
        val = self.valuation()
        iter = ceil(workprec / (val - 1/(parent._p-1))) + 1
        for i in range(1,iter):
            pow *= self / parent(i)
            res += pow
            res = res.truncate(workprec)
        return res

    def teichmuller(self,workprec=Infinity):
        if workprec is Infinity:
            raise ApproximationError("unable compute Teichmuller lift to infinite precision")
        res = self.truncate(1)
        p = self.parent()._p
        for i in range(2,workprec+1):
            res = res ** p
            res = res.truncate(i)
        return res
