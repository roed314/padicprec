from sage.structure.sage_object import SageObject
from sage.misc.classcall_metaclass import ClasscallMetaclass

class PrecisionMode(SageObject):
    __metaclass__ = ClasscallMetaclass
    _mode = [ None for _ in range(3) ]

    @staticmethod
    def __classcall__(cls,mode):
        if isinstance(mode,PrecisionMode):
            return mode
        if not isinstance(mode,str):
            raise TypeError("mode must be an instance of PrecisionMode or a string")
        mode = mode.lower()
        if mode == 'disable': mode = 0
        elif mode == 'individual': mode = 1
        elif mode == 'collective': mode = 2
        else:
            raise ValueError
        if cls._mode[mode] is None:
            cls._mode[mode] = cls.__new__(PrecisionMode)
            cls._mode[mode]._mode = mode
        return cls._mode[mode]

    def _repr_(self):
        if self._mode == 0:
            return "Disable"
        elif self._mode == 1:
            return "Individual"
        elif self._mode == 2:
            return "Collective"

PrecisionMode_Disable = PrecisionMode("disable")
PrecisionMode_Individual = PrecisionMode("individual")
PrecisionMode_Collective = PrecisionMode("collective")


class LazyMode(SageObject):
    __metaclass__ = ClasscallMetaclass
    _mode = [ None for _ in range(4) ]

    @staticmethod
    def __classcall__(cls,mode):
        if isinstance(mode,LazyMode):
            return mode
        if not isinstance(mode,str):
            raise TypeError("mode must be an instance of LazyMode or a string")
        mode = mode.lower()
        if mode == 'disable': mode = 0
        elif mode == 'asap': mode = 1
        elif mode == 'print' or mode == 'printing': mode = 2
        elif mode == 'lazy' or mode == 'never': mode = 3
        else:
            raise ValueError
        if cls._mode[mode] is None:
            cls._mode[mode] = cls.__new__(LazyMode)
            cls._mode[mode]._mode = mode
        return cls._mode[mode]

    def _repr_(self):
        if self._mode == 0:
            return "Disable"
        elif self._mode == 1:
            return "Evaluate as soon as possible"
        elif self._mode == 2:
            return "Evaluate before printing"
        elif self._mode == 3:
            return "Really lazy (evaluate only when necessary)"


LazyMode_Disable = LazyMode("disable")
LazyMode_ASAP = LazyMode("asap")
LazyMode_Printing = LazyMode("printing")
LazyMode_Lazy = LazyMode("lazy")
