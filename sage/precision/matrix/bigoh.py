from sage.rings.infinity import Infinity

from precision.bigoh import BigOh, ExactBigOh, FlatBigOh, JaggedBigOh, LatticeBigOh


class BigOh_matrix(BigOh):
    def __mul__(self,other):
        if self.parent().base_ring() != other.parent().base_ring():
            # should implement coercion here
            raise NotImplementedError
        if self.__class__ is other.__class__:
            return self._mul_(other)

        # Here, we do some convertions
        parent = self.parent()
        path = [ ]
        try:
            # Todo: find something better
            path = parent._models.shortest_path(self.__class__,other.__class__)
            if len(path) > 0:
                for i in range(1,len(path)):
                    map = parent._models.edge_label(path[i-1],path[i])
                    if map == None:
                        self = path[i](parent,self)
                    else:
                        self = map(parent,self)
                return self._mul_(other)
            path = parent._models.shortest_path(other.__class__,self.__class__)
            if len(path) > 0:
                for i in range(1,len(path)):
                    map = parent._models.edge_label(path[i-1],path[i])
                    if map == None:
                        other = path[i](parent,other)
                    else:
                        other = map(parent,other)
                return self._mul_(other)
        except RuntimeError:
            pass
        raise TypeError

    def nrows(self):
        return self.parent()._ambiant_space.nrows()

    def ncols(self):
        return self.parent()._ambiant_space.ncols()

    def _repr_(self,model):
        from sage.functions.other import floor
        l = [ ]
        maxlen = self.ncols() * [0]
        for i in range(self.nrows()):
            lp = [ ]
            for j in range(self.ncols()):
                s = self[i,j].__repr__()
                if len(s) > maxlen[j]: maxlen[j] = len(s)
                lp.append(s)
            l.append(lp)
        s = ""
        for i in range(self.nrows()):
            s += "["
            lp = l[i]
            for j in range(self.ncols()):
                length = len(lp[j]) 
                nbspace = floor((maxlen[j]-length) / 2)
                s += (nbspace + 1)*" " + lp[j] + (maxlen[j] - length - nbspace + 1)*" "
            s += "]\n"
        return s[:-1]


class FlatBigOh_matrix(FlatBigOh, BigOh_matrix):
    def _repr_(self,model):
        return BigOh_matrix._repr_(self,model)


class JaggedBigOh_matrix(JaggedBigOh, BigOh_matrix):
    def _repr_(self,model):
        return BigOh_matrix._repr_(self,model)

    def _mul_(self,other):
        if self.ncols() != other.nrows():
            raise TypeError("incompatible dimensions")
        l = [ ]
        for i in range(self.nrows()):
            for j in range(other.ncols()):
                precision = self._base_exactprec
                for k in range(self.ncols()):
                    precision += self[i,k] * other[k,j]
                l.append(precision)
        from matrix_space import MatrixSpace_inexact
        MSpace = MatrixSpace_inexact(self.parent().base_ring(), self.nrows(), other.ncols())
        parentprec = MSpace.precision()
        return self.__class__(parentprec,l)


class LatticeBigOh_matrix(LatticeBigOh, BigOh_matrix):
    def _repr_(self,model):
        return BigOh_matrix._repr_(self,model)
