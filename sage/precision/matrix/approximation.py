from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.mutability import Mutability
from sage.matrix.matrix_space import MatrixSpace
from sage.matrix.matrix_generic_dense import Matrix_generic_dense

from precision.approximation import Approximation

class MatrixSpaceApprox(MatrixSpace):
    def __init__(self, *args, **kwargs):
        MatrixSpace.__init__(self, *args, **kwargs)
        from lazy import LazyApproximation_matrix
        self._lazy_class = LazyApproximation_matrix
        self._zero = self(0)
        self._lazy_zero = self._lazy_class(self,[])

    def _get_matrix_class(self):
        return MatrixApprox

    def indices_basis(self):
        return [ (i,j) for i in range(self.nrows()) for j in range(self.ncols()) ]


class MatrixApprox(Matrix_generic_dense, Approximation):
    def __init__(self, parent, entries, copy, coerce):
        Matrix_generic_dense.__init__(self, parent, entries, copy, coerce)
        Approximation.__init__(self,parent)
        self._length = parent.dimension()
        self.set_immutable()

    def matrix_space(self, nrows=None, ncols=None, sparse=None):
        if nrows is None:
            nrows = self.nrows()
        if ncols is None:
            ncols = self.ncols()
        if sparse is None:
            sparse = self.is_sparse()
        return MatrixSpaceApprox(self.base_ring(), nrows, ncols, sparse=sparse)

    def truncate(self,workprec):
        return self.__class__(self.parent(), [ c.truncate(workprec) for c in self.list() ], False, False)

    def parenthesis_level(self):
        return 3

    def __repr__(self):
        from printing import repr_matrix
        return repr_matrix([ [ self[i,j] for j in range(self.ncols()) ] for i in range(self.nrows()) ])

    #def valuation_entries(self,lazylimit=None):
    #    return [ self[i,j].valuation() for i in range(self.nrows()) for j in range(self.ncols()) ]

    def _getitem_by_num(self,n):
        ncols = self.ncols()
        j = n % ncols
        i = (n-j) / ncols
        return self[i,j]

    # Grrr
    def _mul_(self,other):
        res = Matrix_generic_dense._mul_(self,other)
        Approximation.__init__(res,res.parent())
        res.set_immutable()
        return res
