from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.graphs.digraph import DiGraph

from precision.exception import PrecisionError
from precision.parent_precision import ParentBigOh
from precision.parent_inexact import Parent_inexact
from bigoh import FlatGroupBigOh, JaggedGroupBigOh, LatticeGroupBigOh


class _ParentGroupBigOh(ParentBigOh,UniqueRepresentation):
    def __init__(self):
        self._models = DiGraph([[FlatGroupBigOh, JaggedGroupBigOh, LatticeGroupBigOh],lambda i,j: False],loops=False)
        self._default_model = LatticeGroupBigOh
        Parent.__init__(self)

    def set_default_model(self,model):
        if self._models.has_vertex(model):
            self._default_model = model
        else:
            raise ValueError

    def default_model(self):
        return self._default_model

    def uniformizer_name(self):
        raise AttributeError("This parent do not have a uniformizer")

    def ambiant_space(self):
        raise AttributeError("This parent do not have an ambiant space")
    
    def __repr__(self):
        return "Parent for all group precisions"


class _ParentGroupElements_inexact(Parent_inexact,UniqueRepresentation):
    def __init__(self):
        Parent.__init__(self)
        self._precision = ParentGroupBigOh

    def precision_mode(self):
        raise PrecisionError("The precision mode is defined individually for each group")

    def set_precision_mode(self,prec):
        raise PrecisionError("The precision mode is defined individually for each group")

    def __repr__(self):
        return "Parent for all groups of inexact elements"


ParentGroupBigOh = _ParentGroupBigOh()
ParentGroupElements_inexact = _ParentGroupElements_inexact()
