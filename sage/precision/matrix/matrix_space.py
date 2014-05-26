from sage.structure.unique_representation import UniqueRepresentation

from sage.rings.infinity import Infinity
from sage.rings.ring import Ring
from sage.graphs.digraph import DiGraph

from sage.rings.morphism import RingHomomorphism_coercion

from precision.parent_inexact import Parent_inexact
from approximation import MatrixSpaceApprox
from bigoh import FlatBigOh_matrix, JaggedBigOh_matrix, LatticeBigOh_matrix


class BaseringInjection(RingHomomorphism_coercion):
    def __init__(self,R,matrix_class):
        RingHomomorphism_coercion.__init__(self, R.base_ring().Hom(R), check=False)
        self._matrix_class = matrix_class

    def _call_(self,x):
        return self._matrix_class(self.codomain(), x)



class MatrixSpace_inexact(Parent_inexact,Ring,UniqueRepresentation):
    @staticmethod
    def __classcall__(cls, base_ring, nrows, ncols=None, element_class=None):
        if element_class is None:
            from matrix import Matrix_inexact
            element_class = Matrix_inexact
        if ncols is None:
            ncols = nrows
        return super(MatrixSpace_inexact,cls).__classcall__(cls, base_ring, nrows, ncols, element_class)

    def __init__(self, base_ring, nrows, ncols, element_class):
        self._matrix_class = element_class
        approximation = MatrixSpaceApprox(base_ring.approximation(), nrows, ncols)
        self._nrows = approximation.nrows()
        self._ncols = approximation.ncols()
        models = DiGraph({FlatBigOh_matrix: [JaggedBigOh_matrix]})
        Parent_inexact.__init__(self, base_ring, approximation, models)
        inject = BaseringInjection(self,self._matrix_class)
        self._populate_coercion_lists_(coerce_list=[inject])
        self._exact_precision = JaggedBigOh_matrix(self._precision, Infinity)

    def _element_constructor_(self, *args, **kwargs):
        return self._matrix_class(self, *args, **kwargs)

    def zero(self):
        return self._matrix_class(self,self.base_ring().zero())
    def zero_element(self):
        return self.zero()
    def zero_matrix(self):
        return self.zero()

    def one(self):
        if self.is_square():
            return self._matrix_class(self,self.base_ring().one())
        else:
            raise TypeError("self must be a space of square matrices")
    def one_element(self):
        return self.one()
    def identity_matrix(self):
        return self.one()

    def _an_element_(self):
        return self.zero()

    def indices_basis(self):
        return [ (i,j) for i in range(self._nrows) for j in range(self._ncols) ]

    def is_commutative(self):
        if self._nrows == 1 and self._ncols == 1:
            return self.base_ring().is_commutative()
        return False

    def nrows(self):
        return self._nrows

    def ncols(self):
        return self._ncols

    def is_square(self):
        return self._nrows == self._ncols

    def _repr_(self):
        return "MatrixSpace of %s by %s matrices over %s" % (self._nrows, self._ncols, self.base_ring())
