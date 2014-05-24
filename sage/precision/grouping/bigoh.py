from precision.bigoh import BigOh
from precision.exception import PrecisionError

class GroupBigOh(BigOh):
    def __init__(self,parent,group):
        self._group = group
        from parents import ParentGroupBigOh
        BigOh.__init__(self, ParentGroupBigOh)

    def __repr__(self):
        return "BigOh for group %s" % self._group

    def __add__(self,other):
        raise PrecisionError("Can't add two GroupBigOh")

    def add_element(self,key,prec):
        #print "Add element %s" % key()
        pass

    def delete_element(self,key):
        #print "Delete element %s" % key()
        pass

    def precision(self,key):
        # key can be a single element or a list of keys
        raise NotImplementedError

    def to_workprec(self):
        return 0


class FlatGroupBigOh(GroupBigOh):
    pass

class JaggedGroupBigOh(GroupBigOh):
    pass

class LatticeGroupBigOh(GroupBigOh):
    pass
