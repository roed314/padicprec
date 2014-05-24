import _weakref as weakref
import inspect
from types import FunctionType
from sage.misc.cachefunc import cached_function

from sage.structure.sage_object import SageObject
from sage.structure.element import Element
from sage.rings.integer import Integer
from sage.rings.infinity import Infinity
from approximation import Approximation
from limits import maxworkprec, maxworkprec_cache

from exception import PrecisionError


class LazyApproximation(Approximation):
    def __init__(self, parent, x, max_workprec=None, **kwargs):
        Approximation.__init__(self, parent)
        self._length = -1
        self._current_workprec = None
        self._current_approximation = None
        self._dict_max_workprec = { }
        if isinstance(x,Approximation) and x.parent() is parent:
            if isinstance(x,LazyApproximation):
                raise TypeError("can't lazy something already lazy")
            self._current_workprec = Infinity
            self._current_approximation = x
            self.repr = lambda: "%s" % x
            self._length = length = x._length
            base = parent.base_ring()
            if length >= 0:
                zero = base._lazy_zero
                lazy_class = base._lazy_class
                @cached_function
                def getitem_by_num(i):
                    if i < length:
                        return lazy_class(base, x._getitem_by_num(i)).new_wrapper()
                    else:
                        return zero
                self._getitem_by_num = getitem_by_num
        elif isinstance(x,list):
            self.name = 'element_from_list'
            if parent.dimension() < Infinity:
                zeros = [ parent.base_ring()._zero ] * (parent.dimension() - len(x))
                self.function = lambda workprec: parent([ item.at_workprec(workprec) for item in x ] + zeros)
            else:
                self.function = lambda workprec: parent([ item.at_workprec(workprec) for item in x ])
            x_wrapper = [ c.new_wrapper() for c in x ]
            self._length = length = len(x)
            zero = parent.base_ring()._lazy_zero
            def getitem_by_num(i):
                if i < length:
                    return x_wrapper[i]
                else:
                    return zero
            self._getitem_by_num = getitem_by_num
            self.args = self._dependances = x_wrapper
            if len(x) == 0:
                if max_workprec is None:
                    max_workprec = 0
            else:
                if max_workprec is None:
                    max_workprec = max([ item.max_workprec() for item in x ])
        elif isinstance(x,dict):
            raise NotImplementedError
        elif isinstance(x, FunctionType):
            self.function = x
            import re
            self.name = re.sub("^(lazy)?_*", "", x.__name__)
            self.name = re.sub("_*$", "", self.name)
        else:
            raise TypeError("x must be an approximation, a list, a dictionary or a function")

        if max_workprec is None:
            self._max_workprec = Infinity
        else:
            self._max_workprec = max_workprec

        if kwargs.has_key('parenthesis_level'):
            self._parenthesis_level = kwargs['parenthesis_level']
            del kwargs['parenthesis_level']
        else:
            self._parenthesis_level = 3
            try:
                if self.name == 'add' or self.name == 'sub' or self.name == 'neg': self._parenthesis_level = 0
                elif self.name == 'mul' or self.name == 'div' or self.name == 'invert': self._parenthesis_level = 1
                elif self.name == 'pow': self._parenthesis_level = 2
            except AttributeError:
                pass

        for key in kwargs:
            if key[:1] != '_':
                self.__dict__[key] = kwargs[key]

    def __hash__(self):
        return id(self)

    def new_wrapper(self,wrapper=None):
        if wrapper is None:
            workprec = self._max_workprec
        else:
            workprec = self._dict_max_workprec[wrapper]
        return LazyApproximationWrapper(self,workprec)

    def __repr__(self, *args, **kwargs):
        return self.repr(*args, **kwargs)

    def _del_max_workprec(self,key):
        d = self._dict_max_workprec
        del d[key]
        if len(d) > 0:
            self._max_workprec = max(d.values())
        self.replace_by_value()

    def max_workprec(self,wrapper=None):
        if wrapper is None:
            return self._max_workprec
        else:
            return self._dict_max_workprec[wrapper]

    def update_max_workprec(self,workprec,wrapper=None):
        d = self._dict_max_workprec
        try:
            d[wrapper] = min(d[wrapper], workprec)
        except KeyError:
            d[wrapper] = workprec
        self._max_workprec = max(d.values())
        try:
            for c in self._dependances:
                c.update_max_workprec(workprec)
        except AttributeError:
            pass

    def current_workprec(self):
        return self._current_workprec

    def parenthesis_level(self):
        return self._parenthesis_level

    def repr(self):
        from printing import repr_parenthesis, repr_operator, repr_exponant, repr_index, repr_function
        name = self.name
        if self.__dict__.has_key('args'):
            args = [ ]
            parenthesis = [ ]
            for arg in self.args:
                level = 0
                if hasattr(arg, 'parenthesis_level'):
                    level = arg.parenthesis_level()
                elif hasattr(arg, 'is_atomic') and self._parenthesis_level < Infinity:
                    if arg.is_atomic(): level = 3
                elif isinstance(arg, Integer):
                    level = 3
                args.append("%s" % arg)
                parenthesis.append(level)
            s = ""
            if name == 'add':
                s = repr_operator(" + ", args)
            elif name == 'sub' and len(args) == 2:
                if parenthesis[1] < 1:
                    args[1] = repr_parenthesis(args[1])
                s = repr_operator(" - ", args)
            elif name == 'neg' and len(args) == 1:
                if parenthesis[0] < 1:
                    args[0] = repr_parenthesis("%s" % args[0])
                s = repr_operator("-", ["", args[0]])
            elif name == 'mul':
                for i in range(len(args)):
                    if parenthesis[i] < 1:
                        args[i] = repr_parenthesis(args[i])
                s = repr_operator("*", args)
            elif name == 'div' and len(args) == 2:
                if parenthesis[0] < 1:
                    args[0] = repr_parenthesis(args[0])
                if parenthesis[1] < 2:
                    args[1] = repr_parenthesis(args[1])
                s = repr_operator("/", args)
            elif name == 'invert' and len(args) == 1:
                if parenthesis[0] < 2:
                    args[0] = repr_parenthesis(args[0])
                s = repr_operator("/", ["1", args[0]])
            elif name == 'pow' and len(args) == 2:
                if parenthesis[0] < 3:
                    args[0] = repr_parenthesis(args[0])
                if parenthesis[1] < 3:
                    args[1] = repr_parenthesis(args[1])
                s = repr_exponant(args[0], args[1])
            elif name == 'getitem' and len(args) == 2:
                if parenthesis[0] < 2:
                    args[0] = repr_parenthesis("%s" % args[0])
                s = repr_index(args[0], args[1])
            else:
                s = repr_function(name, args)
            return s
        else:
            return "%s(...)" % name

    def _getitem_by_num(self,i):
        def lazy_getitem(workprec):
            return self.at_workprec(workprec)._getitem_by_num(i)
        base = self.base_ring()
        indices = self.parent().indices_basis()
        if indices is not None:
            index = indices[i]
        else:
            index = i
        return base._lazy_class(base, lazy_getitem, args=[self,str(index)])

    def parenthesis_level(self):
        return self._parenthesis_level

    def at_workprec(self,workprec=None,wrapper=None):
        if wrapper is None:
            max_workprec = self._max_workprec
        else:
            max_workprec = self._dict_max_workprec[wrapper]
        if workprec is None or workprec > max_workprec:
            workprec = max_workprec
        #if workprec > max_workprec and workprec is not Infinity:
        #    raise PrecisionError("working precision overflow")
        if self._current_workprec is None or workprec > self._current_workprec:
            if workprec > maxworkprec_cache and workprec is not Infinity:
                res = self.function(workprec)
            else:
                self._current_approximation = res = self.function(workprec).truncate(workprec)
                self._current_workprec = workprec
        else:
            res = self._current_approximation.truncate(workprec)
        self.replace_by_value()
        return res

    def replace_by_value(self):
        if self._current_workprec is not None and self._current_workprec >= self._max_workprec:
            unlazy = self._current_approximation
            self.repr = lambda: unlazy.__repr__()
            self.function = None
            self.args = None
            self._unlazy = lambda: None
            if self._length < 0 and unlazy._length >= 0:
                length = self._length = unlazy._length
                parent = self.parent().base_ring()
                zero = parent._lazy_zero
                lazy_class = parent._lazy_class
                @cached_function
                def getitem_by_num(i):
                    if i < length:
                        return lazy_class(parent, unlazy._getitem_by_num(i)).new_wrapper()
                    else:
                        return zero
                self._getitem_by_num = getitem_by_num
            self._parenthesis_level = unlazy.parenthesis_level()

    def is_lazy(self):
        return True

    def length(self,wrapper=None):
        if wrapper is None:
            max_workprec = self._max_workprec
        else:
            max_workprec = self._dict_max_workprec[wrapper]
        if self._current_workprec >= max_workprec:
            return self._current_approximation._length
        else:
            return self._length

    def is_evaluated(self,wrapper=None):
        if wrapper is None:
            max_workprec = self._max_workprec
        else:
            max_workprec = self._dict_max_workprec[wrapper]
        return self._current_workprec >= max_workprec

    def is_evaluated_zero(self,wrapper=None):
        if wrapper is None:
            max_workprec = self._max_workprec
        else:
            max_workprec = self._dict_max_workprec[wrapper]
        if self._current_workprec < max_workprec or self._current_approximation is None:
            return False
        return self._current_approximation.truncate(max_workprec).is_zero()

    # Usual operations

    def _add_(self,other,expected_valuation=0):
        if self.is_evaluated_zero():
            return other
        if other.is_evaluated_zero():
            return self
        parent = self.parent()
        length = max(self._length, other._length)
        if self._length >= 0 and other._length >= 0:
            res = self.__class__(parent, [ self._getitem_by_num(i)._add_(other._getitem_by_num(i), expected_valuation) for i in range(length) ])
        else:
            def lazy_add(workprec):
                return self.at_workprec(workprec)._add_(other.at_workprec(workprec))
            res = self.__class__(parent, lazy_add, args=[self,other], expected_valuation=expected_valuation)
            if self._length >= 0 or other._length >= 0:
                parentbase = parent.base_ring()
                classbase = parentbase._lazy_class
                @cached_function
                def getitem_by_num(i):
                    itemself = self._getitem_by_num(i)
                    itemother = other._getitem_by_num(i)
                    def lazy_add(workprec):
                        return itemself.at_workprec(workprec)._add_(itemother.at_workprec(workprec))
                    return classbase(parentbase, lazy_add, args=[itemself, itemother])
                res._getitem_by_num = getitem_by_num
        return res

    def _sub_(self,other,expected_valuation=0):
        if other.is_evaluated_zero():
            return self
        if self.is_evaluated_zero():
            return other.__neg__(expected_valuation=expected_valuation)
        parent = self.parent()
        length = max(self._length, other._length)
        if self._length >= 0 or other._length >= 0:
            # should we modify lazy_options?
            res = self.__class__(parent, [ self._getitem_by_num(i)._sub_(other._getitem_by_num(i), expected_valuation=expected_valuation) for i in range(length) ])
        else:
            def lazy_sub(workprec):
                return self.at_workprec(workprec)._sub_(other.at_workprec(workprec))
            res = self.__class__(parent, lazy_sub, args=[self,other], expected_valuation=expected_valuation)
            if self._length >= 0 or other._length >= 0:
                parentbase = parent.base_ring()
                classbase = parentbase._lazy_class
                @cached_function
                def getitem_by_num(i):
                    itemself = self._getitem_by_num(i)
                    itemother = other._getitem_by_num(i)
                    def lazy_sub(workprec):
                        return itemself.at_workprec(workprec)._sub_(itemother.at_workprec(workprec))
                    # What should we do with lazy_options here?
                    return classbase(parentbase, lazy_sub, args=[itemself, itemother])
                res._getitem_by_num = getitem_by_num
        return res

    def __neg__(self,expected_valuation=0):
        parent = self.parent()
        if self.is_evaluated_zero():
            return parent._lazy_zero
        length = self._length
        if self._length >= 0:
            # should we modify lazy_options?
            res = self.__class__(parent, [ self._getitem_by_num(i).__neg__(expected_valuation=expected_valuation) for i in range(length) ])
        else:
            def lazy_neg(workprec):
                return self.at_workprec(workprec).__neg__()
            res = self.__class__(parent, lazy_neg, args=[self], expected_valuation=expected_valuation)
        return res

    def _rmul_(self,coeff,workprec_shifts,expected_valuation=0):
        parent = self.parent()
        if self.is_evaluated_zero() or coeff.is_evaluated_zero():
            return parent._lazy_zero
        length = self._length
        if self._length >= 0:
            # should we modify lazy_options?
            res = self.__class__(parent, [ self._getitem_by_num(i)._rmul_(coeff, workprec_shifts, expected_valuation=expected_valuation) for i in range(length) ])
        else:
            (shiftself, shiftcoeff) = workprec_shifts
            def lazy_mul(workprec):
                return self.at_workprec(workprec + shiftself)._rmul_(coeff.at_workprec(workprec + shiftcoeff))
            res = self.__class__(parent, lazy_mul, args=[coeff,self], expected_valuation=expected_valuation)
        return res

    def _lmul_(self,coeff,workprec_shifts,expected_valuation=0):
        parent = self.parent()
        if self.is_evaluated_zero() or coeff.is_evaluated_zero():
            return parent._lazy_zero
        length = self._length
        if self._length >= 0:
            # should we modify lazy_options?
            res = self.__class__(parent, [ self._getitem_by_num(i)._lmul_(coeff, workprec_shifts, expected_valuation=expected_valuation) for i in range(length) ])
        else:
            (shiftself, shiftcoeff) = workprec_shifts
            def lazy_mul(workprec):
                return self.at_workprec(workprec + shiftself)._lmul_(coeff.at_workprec(workprec + shiftcoeff))
            res = self.__class__(parent, lazy_mul, args=[self,coeff], expected_valuation=expected_valuation)
        return res

    def _mul_(self,other,workprec_shifts,expected_valuation=0):
        parent = self.parent()
        if self.is_evaluated_zero() or other.is_evaluated_zero():
            return parent._lazy_zero
        (shiftself, shiftother) = workprec_shifts
        def lazy_mul(workprec):
            return self.at_workprec(workprec + shiftself)._mul_(other.at_workprec(workprec + shiftother))
        return self.__class__(parent, lazy_mul, args=[self,other], expected_valuation=expected_valuation)

    def _pow_(self,exp,workprec_transformation,expected_valuation=0):
        parent = self.parent()
        def lazy_pow(workprec):
            return self.at_workprec(workprec_transformation(workprec)).__pow__(exp)
        return self.__class__(parent, lazy_pow, args=[self,exp], expected_valuation=expected_valuation)


wrappers = { }
def delete_max_workprec(ref):
    wrappers[ref]._del_max_workprec(ref)
    del wrappers[ref]

class LazyApproximationWrapper(LazyApproximation):
    def __init__(self,instance,workprec=None):
        if isinstance(instance, LazyApproximationWrapper):
            raise TypeError
        Element.__init__(self, instance.parent())
        self._instance = instance
        ref = weakref.ref(self, delete_max_workprec)
        wrappers[ref] = instance
        self._ref = ref
        if workprec is not None:
            instance.update_max_workprec(workprec,wrapper=ref)

    def __hash__(self):
        return id(self)

    def __repr__(self):
        return self.repr()

    def __getattribute__(self, name):
        instance = LazyApproximation.__getattribute__(self, '_instance')
        attr = instance.__getattribute__(name)
        try:
            vars = inspect.getargspec(attr)[0]
            vars.index('wrapper')
        except (TypeError,ValueError):
            return attr
        ref = Approximation.__getattribute__(self, '_ref')
        def modified_attr(*args, **kwargs):
            kwargs['wrapper'] = ref
            return attr(*args, **kwargs)
        return modified_attr

    def __getitem__(self,i):
        return self.__getitem__(i)

    # for Tab-completion
    def __dir__(self):
        return dir(LazyApproximation.__getattribute__(self, '_instance'))


def lazy_method(method, lazy_name=None):
    if lazy_name is None:
        lazy_name = method.__name__
    def lazy_method(*lazy_args, **options):
        if options.has_key('lazy_options'):
            lazy_options = options['lazy_options']
            del options['lazy_options']
        else:
            lazy_options = { }
        nbargs = len(lazy_args)
        if lazy_options.has_key('workprec_transformations'):
            workprec_transformations = lazy_options['workprec_transformations']
            del lazy_options['workprec_transformations']
        else:
            workprec_transformations = nbargs * [ None ]
        def function(workprec):
            args = [ ]
            for i in range(nbargs):
                if isinstance(lazy_args[i], LazyApproximation):
                    if workprec_transformations[i] is not None:
                        wprec = workprec_transformations[i](workprec)
                    else:
                        wprec = workprec
                    args.append(lazy_args[i].at_workprec(wprec))
                else:
                    args.append(lazy_args[i])
            opt = { }
            for key in options:
                value = options[key]
                if isinstance(value, FunctionType):
                    opt[key] = value(workprec)
                else:
                    opt[key] = value
            return method(*args, workprec=workprec, **opt)
        function.__name__ = lazy_name
        if lazy_options.has_key('parent_answer'):
            parent = lazy_options['parent_answer']
            del lazy_options['parent_answer']
        else:
            parent = lazy_args[0].parent()
        lazy_options['args'] = lazy_args
        return parent._lazy_class(parent, function, **lazy_options)
    lazy_method.__name__ = method.__name__
    return lazy_method
