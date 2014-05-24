import _weakref as weakref
from sage.structure.sage_object import SageObject
from sage.rings.infinity import Infinity
from sage.rings.integer_ring import ZZ

from precision.element_inexact import Element_inexact

from mode import PrecisionMode, LazyMode
from mode import PrecisionMode_Collective

from exception import GroupError

class GroupElements_inexact(SageObject):
    def __init__(self, label=None, precision_mode=None, capped_relative=None, lazy_mode=None, remember_members=None):
        if label is not None and not isinstance(label,str):
            raise TypeError("label must be None or a string")
        self._label = label
        try:
            from context import stack
            current_group = stack[-1]
        except ImportError:
            pass
        if precision_mode is None:
            self._precision_mode = current_group.precision_mode()
        else:
            self._precision_mode = PrecisionMode(precision_mode)
        if capped_relative is None:
            self._capped_relative = current_group.capped_relative()
        elif capped_relative is Infinity:
            self._capped_relative = Infinity
        else:
            self._capped_relative = ZZ(capped_relative)
        if lazy_mode is None:
            self._lazy_mode = current_group.lazy_mode()
        else:
            self._lazy_mode = LazyMode(lazy_mode)
        if remember_members is None:
            self._remember = self._precision_mode is PrecisionMode_Collective
        else:
            self._remember = remember_members
        if self._remember:
            self._members = set()

    def label(self):
        return self.label

    def precision_mode(self):
        return self._precision_mode

    def capped_relative(self):
        return self._capped_relative

    def lazy_mode(self):
        return self._lazy_mode

    def remember_members(self):
        return self._remember

    def __repr__(self):
        import re
        if self._label is not None:
            return "Group '%s'" % self._label
        if self._remember:
            members = self.members()
            if len(members) > 0:
                s = "Group consisting of the following element(s):"
                for element in self.members():
                    s += "\n * " + re.sub("\n", "\n   ", element.__repr__())
            else:
                s = "Empty group of inexact elements"
            return s
        else:
            return "Group of inexact elements"

    def members(self,safemode=False):
        if not self._remember:
            if safemode:
                return
            else:
                raise GroupError("This group does not remember its members")
        return [ ref() for ref in self._members ]

    def register_member(self,member,safemode=False):
        if not self._remember:
            if safemode:
                return
            else:
                raise GroupError("This group does not remember its members")
        if not isinstance(member, Element_inexact):
            raise TypeError("member must be an instance of Element_inexact")
        ref = weakref.ref(member, self.__delete_member_by_ref)
        self._members.add(ref)
        if self._precision_mode is PrecisionMode_Collective:
            # do something here
            pass

    def discard_member(self,member,safemode=False):
        if not self._remember:
            if safemode:
                return
            else:
                raise GroupError("This group does not remember its members")
        self.__delete_member_by_ref(weakref.ref(member))

    def __delete_member_by_ref(self,ref):
        self._members.remove(ref)
        if self._precision_mode is PrecisionMode_Collective: 
            # do something here
            pass
