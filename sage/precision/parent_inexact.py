from sage.structure.unique_representation import UniqueRepresentation

from sage.structure.parent import Parent
from sage.rings.integer_ring import ZZ
from sage.rings.infinity import Infinity

from parent_precision import ParentBigOh
#from precision_mode import *
from bigoh import ExactBigOh, ValuationBigOh, FlatBigOh, FlatIntervalBigOh, JaggedBigOh
from element_inexact import Element_inexact

from exception import PrecisionError


class Parent_inexact(Parent,UniqueRepresentation):
    def __init__(self, base, parent_approximation, models=None, exact=Infinity,
                 category=None, element_constructor=None, gens=None, names=None, normalize=True, facade=None):
        Parent.__init__(self,base,category=category,element_constructor=element_constructor,
                        gens=gens, names=names, normalize=normalize, facade=facade)
        self._approximation = parent_approximation
        self._precision = ParentBigOh(self)
        if models is None:
            dimension = self.dimension()
            if dimension == 1:
                self._precision.add_model(ValuationBigOh)
            else:
                if self.dimension() < Infinity:
                    self._precision.add_model(FlatBigOh)
                self._precision.add_model(JaggedBigOh)
                try:
                    if self.indices_basis() is None:
                        self._precision.add_model(FlatIntervalBigOh)
                        self._precision.add_modelconversion(FlatIntervalBigOh,JaggedBigOh)
                except NotImplementedError:
                    pass
                self._precision.add_modelconversion(FlatBigOh,JaggedBigOh,safemode=True)
                self._precision.add_modelconversion(FlatBigOh,FlatIntervalBigOh,safemode=True)
        else:
            for model in models.vertices():
                self._precision.add_model(model)
            for (mfrom,mto,map) in models.edges():
                self._precision.add_modelconversion(mfrom,mto,map)
        try:
            self._exact_precision = self._precision(exact)
        except PrecisionError:
            self._exact_precision = ExactBigOh(self._precision)

    def __hash__(self):
        return id(self)

    def uniformizer_name(self):
        base = self._base
        if hasattr(base, "uniformizer_name"):
            return base.uniformizer_name()
        elif hasattr(base, "uniformizer"):
            import re
            return re.sub("\s*\+\s*O(.*)", "", base.uniformizer().__repr__())
        else:
            raise NotImplementedError

    def _repr_(self):
        s = self._approximation.__repr__()
        i = s.rfind("over ")
        if i == -1: s = ""
        else: s = s[:i+5]
        s += self._base.__repr__()
        return s

    def approximation(self):
        return self._approximation

    def default_working_precision(self):
        return self._default_workprec

    def is_exact(self):
        return False

    def random_element(self,*args,**kwargs):
        dimension = self.dimension()
        if dimension is Infinity:
            raise NotImplementedError
        base = self.base_ring()
        coefficients = [ base.random_element(*args,**kwargs) for _ in range(dimension) ]
        return self(coefficients)

    def dimension(self):
        try:
            return self._approximation.dimension()
        except AttributeError:
            return Infinity

    def indices_basis(self):
        """
        None: the basis is indexed by usual natural integers (including 0, as always in Python)
        """
        pass

    def precision_models(self,graph=False):
        return self._precision.models(graph)

    def precision(self):
        return self._precision
