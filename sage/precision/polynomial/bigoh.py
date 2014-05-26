from sage.rings.integer import Integer
from sage.rings.infinity import Infinity

from polynomial_element import Polynomial_inexact
from precision.bigoh import BigOh, FlatIntervalBigOh, JaggedBigOh, LatticeBigOh
from precision.exception import PrecisionError


class BigOh_polynomial(BigOh):
    def _repr_(self,model):
        s = ""
        m = self.last()
        name = self.parent().ambiant_space().variable_name()
        for n in range(m-1,-1,-1):
            prec = self[n]
            if prec.is_exact():
                continue
            s += " + %s" % prec
            if n > 0:
                s += "*%s" % name
            if n > 1:
                s += "^%s" % n
        if s == "": return "Exact"
        return s[3:]


class FlatIntervalBigOh_polynomial(FlatIntervalBigOh, BigOh_polynomial):
    def __init__(self, *args, **kwargs):
        FlatIntervalBigOh.__init__(self, *args, **kwargs)
        if self.last() is Infinity:
            raise PrecisionError("Infinite precision for polynomials not allowed")

    def _repr_(self,model):
        return BigOh_polynomial._repr_(self,model)


class JaggedBigOh_polynomial(JaggedBigOh, BigOh_polynomial):
    def _mul_(self,other):
        precision = [ self._base_exactprec ] * (other.last() + self.last() - 1)
        for i in range(self.last()):
            for j in range(other.last()):
                precision[i+j] += self[i] * other[j]
        return self.__class__(self.parent(), precision)

    def _repr_(self,model):
        return BigOh_polynomial._repr_(self,model)


class NewtonBigOh_polynomial(BigOh_polynomial):
    def _init_(self,parent,prec,check=True):
        from sage.geometry.newton_polygon import NewtonPolygon, NewtonPolygon_element
        if isinstance(prec,NewtonPolygon_element):
            self._polygon = prec
        else:
            if not isinstance(prec,list):
                prec = [ prec ]
                #raise TypeError("prec must be either a Newton polygon, a list of coordinates (x,y) or a list of precisions")
            dimension = parent.dimension()
            vertices = [ ]
            if len(prec) > 0:
                if isinstance(prec[0],tuple):
                    for p in prec:
                        if isinstance(p,tuple) and (len(p) == 2):
                            if (p[0] < 0) or (p[0] >= dimension):
                                raise IndexError
                            val = p[1]
                            if isinstance(val,Integer) or val is Infinity:
                                vertices.append((p[0],val))
                            else:
                                try:
                                    val = self._baseprec(val).flat_below()
                                except TypeError:
                                    raise TypeError("prec must be either a Newton polygon, a list of coordinates (x,y) or a list of precisions")
                                vertices.append((p[0],val))
                        else:
                            raise TypeError("prec must be either a Newton polygon, a list of coordinates (x,y) or a list of precisions")
                else:
                    if len(prec) > dimension:
                        raise ValueError("list too long")
                    for i in range(len(prec)):
                        val = prec[i]
                        if isinstance(val,Integer) or val is Infinity:
                            vertices.append((i,val))
                        else:
                            try:
                                val = self._baseprec(val).flat_below()
                            except TypeError:
                                raise TypeError("prec must be either a Newton polygon, a list of coordinates (x,y) or a list of precisions")
                            vertices.append((i,val))
            self._polygon = NewtonPolygon(vertices)

    def _convert_(self,parent,prec):
        vertices = [ ]
        for i in range(prec.last()):
            vertices.append((i,prec[i].flat_below()))
        from sage.geometry.newton_polygon import NewtonPolygon
        self._polygon = NewtonPolygon(vertices)

    def is_exact(self):
        return len(self._polygon.vertices()) == 0

    def _getitem_by_num(self,i):
        from sage.functions.other import ceil
        val = self._polygon(i)
        if val is Infinity:
            return self._base_exactprec
        else:
            return self._baseprec(ceil(val))

    def _add_(self,other):
        return self.__class__(self.parent(), self._polygon + other._polygon)

    def _mul_(self,other):
        return self.__class__(self.parent(), self._polygon * other._polygon)

    def _lmul_(self,coeff):
        return self.__class__(self.parent(), self._polygon.__lshift__(coeff.flat_below()))

    def _rmul_(self,coeff):
        return self.__class__(self.parent(), self._polygon.__lshift__(coeff.flat_below()))

    def __pow__(self,exp,ignored=None):
        return self.__class__(self.parent(), self._polygon ** exp)

    def first(self):
        vertices = self._polygon.vertices()
        if len(vertices) == 0:
            return self.parent().dimension()
        else:
            return vertices[0][0]

    def last(self):
        vertices = self._polygon.vertices()
        if len(vertices) == 0:
            return 0
        else:
            return vertices[-1][0] + 1

    def __lshift__(self,i):
        return self.__class__(self.parent(), self._polygon.__lshift__(i))

    def __rshift__(self,i):
        return self.__class__(self.parent(), self._polygon.__rshift__(i))

    def flat_below(self):
        vertices = self._polygon.vertices()
        if len(vertices) == 0:
            return Infinity
        else:
            return min([ v[1] for v in vertices ])

    def to_workprec(self,last=None):
        vertices = self._polygon.vertices()
        if len(vertices) == 0:
            return Infinity
        last_x = vertices[-1][0] + 1
        if last is None:
            last = last_x
        if last > len(vertices):
            return Infinity
        elif last < len(vertices):
            raise ValueError("still inexact BigOh's after last")  # is it really an error?
        else:
            return max([ v[1] for v in vertices ])

    def equals(self,other):
        return self._polygon == other._polygon

    def contains(self,other):
        return self._polygon <= other._polygon

    def _repr_(self,model):
        return BigOh_polynomial._repr_(self,model)


class LatticeBigOh_polynomial(LatticeBigOh, BigOh_polynomial):
    def _repr_(self,model):
        return BigOh_polynomial._repr_(self,model)
