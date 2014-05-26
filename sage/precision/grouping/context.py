from sage.structure.sage_object import SageObject 
from sage.rings.infinity import Infinity

from mode import PrecisionMode, LazyMode 
from group import GroupElements_inexact 
from precision.element_inexact import Element_inexact

from exception import GroupError

default_group = GroupElements_inexact('Default', 'individual', Infinity, 'asap', False) 
stack = [ default_group ]


def set_default_context(label=None, precision_mode=None, capped_relative=None, lazy_mode=None, remember_members=None):
    if label is None:
        label = 'Default'
    global stack
    if len(stack) > 1:
        raise GroupError
    stack = [ GroupElements_inexact(label, precision_mode, capped_relative, lazy_mode, remember_members) ]

def unset_default_context():
    global stack
    if len(stack) > 1:
        raise GroupError
    stack = [ default_group ]


class Context_general(object):
    def __init__(self, label=None, precision_mode=None, capped_relative=None, lazy_mode=None):
        self._precision_mode = precision_mode
        self._capped_relative = capped_relative
        self._lazy_mode = lazy_mode
        if label is None:
            s = "Context: "
            if precision_mode is not None:
                s += "precision_mode = %s, " % self._precision_mode
            if capped_relative is not None:
                s += "capped_relative = %s, " % self._capped_relative
            if lazy_mode is not None:
                s += "lazy_mode = %s, " % self._lazy_mode
            if s[-2:] == ', ':
                s = s[:-2]
            self._label = s
        else:
            self._label = label
        self._lengths = [ ]
        self._count = 0

    def new_group(self):
        label = self._label
        self._count += 1
        if self._count > 1:
            label += " - %s" % self._count
        return GroupElements_inexact(label, self._precision_mode, self._capped_relative, self._lazy_mode)

    def push_to_stack(self):
        group = self.new_group()
        stack.append(group)
        self._lengths.append(len(stack))
        return group

    def pop_stack(self):
        if len(stack) != self._lengths[-1]:
            raise GroupError
        self._lengths.pop()
        return stack.pop()

    def __repr__(self):
        return label


class Context(Context_general):
    def __enter__(self):
        return self.push_to_stack()

    def __exit__(self,type,value,traceback):
        self.pop_stack()
        return False


class ContextDecorator(Context_general):
    def __init__(self, function, **kwargs):
        self._function = function
        Context_general.__init__(self, **kwargs)

    def __call__(self, *args, **kwargs):
        self.push_to_stack()
        modified_args = self.deep_change_group(args)
        modified_kwargs = self.deep_change_group(kwargs)
        value = self._function(*modified_args, **modified_kwargs)
        self.pop_stack()
        return self.deep_change_group(value)

    def deep_change_group(self,value,group=None):
        if isinstance(value,Element_inexact):
            return value.change_group(group)
        # Probably, there is something better
        elif isinstance(value,dict):
            return { key: self.deep_change_group(value[key],group) for key in value }
        elif isinstance(value,list):
            return [ self.deep_change_group(v,group) for v in value ]
        elif isinstance(value,tuple):
            return tuple([ self.deep_change_group(v,group) for v in value ])
        else:
            return value

    def __repr__(self):
        s = "%s with decoration %s" % (self._function, self._label)


def context(**kwargs):
    def decorator(function):
        return ContextDecorator(function, **kwargs)
    return decorator
