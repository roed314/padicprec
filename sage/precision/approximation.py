import _weakref as weakref

from sage.structure.element import Element
from sage.rings.infinity import Infinity


class Approximation(Element):
    def __init__(self,parent):
        Element.__init__(self,parent)
        self._length = -1
        self._max_workprec = Infinity

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

    def length(self):
        return self._length

    # for lazy purpose
    def current_workprec(self):
        return Infinity

    def max_workprec(self):
        return self._max_workprec

    def update_max_workprec(self,workprec):
        self._max_workprec = min(self._max_workprec,workprec)

    def at_workprec(self,workprec=None):
        if workprec is None:
            workprec = self._max_workprec
        if workprec is Infinity:
            return self
        else:
            return self.truncate(workprec)

    def parenthesis_level(self):
        raise NotImplementedError

    def truncate(self,workprec):
        raise NotImplementedError

    def is_lazy(self):
        return False

    def is_evaluated_zero(self):
        return self.is_zero()

    def is_evaluated(self):
        return True
