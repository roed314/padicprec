from sage.rings.ring import Ring
from sage.graphs.digraph import DiGraph
from sage.rings.morphism import RingHomomorphism_coercion

from precision.parent_inexact import Parent_inexact
from approximation import PolynomialRingApprox
from bigoh import NewtonBigOh_polynomial, FlatIntervalBigOh_polynomial, JaggedBigOh_polynomial, LatticeBigOh_polynomial


class BaseringInjection(RingHomomorphism_coercion):
    def __init__(self,R,polynomial_class):
        RingHomomorphism_coercion.__init__(self, R.base_ring().Hom(R), check=False)
        self._polynomial_class = polynomial_class

    def _call_(self,x):
        return self._polynomial_class(self.codomain(), [x])

    def _call_with_args(self,x,args,kwargs):
        if len(args) > 0:
            prec = args[0]
            args = list(args[1:])
            return self._polynomial_class(self.codomain(), [x], [prec], *args, **kwargs)
        else:
            return self._polynomial_class(self.codomain(), [x], **kwargs)


class PolynomialRing_inexact(Parent_inexact,Ring):
    def __init__(self, base_ring, name=None, names=None, element_class=None):
        if names is None:
            names = name
        if names is None:
            raise ValueError("You must specify a variable name")
        if element_class:
            self._polynomial_class = element_class
        else:
            from polynomial_element import Polynomial_inexact
            self._polynomial_class = Polynomial_inexact
        approximation = PolynomialRingApprox(base_ring.approximation(), name=names)
        models = DiGraph({FlatIntervalBigOh_polynomial: [NewtonBigOh_polynomial], JaggedBigOh_polynomial: [NewtonBigOh_polynomial]})
        Parent_inexact.__init__(self, base_ring, approximation, models, exact=[], names=names)
        self._generator = self._polynomial_class(self,[base_ring.zero(),base_ring.one()])
        inject = BaseringInjection(self,self._polynomial_class)
        self._populate_coercion_lists_(coerce_list=[inject])

    def _element_constructor_(self,x,*args,**kwargs):
        if isinstance(x,int):
            return self._polynomial_class(self,[self.base_ring()(x)])
        return self._polynomial_class(self,x,*args,**kwargs)

    def variable_name(self):
        return self._names[0]

    def _repr_(self):
        return "Polynomial Ring in %s over %s" % (self.variable_name(), self.base_ring())

    def is_commutative(self):
        return True

    def gen(self):
        return self._generator

    def gens(self):
        return (self._generator,)

    def random_element(self,degree=2,monic=False,*args,**kwargs):
        base = self.base_ring()
        if isinstance(degree, (list, tuple)):
            if len(degree) != 2:
                raise ValueError("degree argument must be an integer or a tuple of 2 integers (min_degree, max_degree)")
            if degree[0] > degree[1]:
                raise ValueError("minimum degree must be less or equal than maximum degree")
            degree = randint(*degree)
        if monic:
            return self([ base.random_element(*args, **kwargs) for _ in range(degree) ] + [ base(1) ])
        else:
            return self([ base.random_element(*args, **kwargs) for _ in range(degree+1) ])
