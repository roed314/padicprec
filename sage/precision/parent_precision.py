from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent

from sage.graphs.digraph import DiGraph

from exception import PrecisionError


class ParentBigOh(Parent,UniqueRepresentation):
    def __init__(self,ambiant):
        try:
            self._uniformizer = ambiant.uniformizer_name()
        except NotImplementedError:
            raise TypeError("Impossible to determine the name of the uniformizer")
        self._ambiant_space = ambiant
        self._models = DiGraph(loops=False)
        self._default_model = None
        Parent.__init__(self,ambiant.base_ring())
        if self.base_ring() is None:
            self._base_precision = None
        else:
            if self.base_ring() == ambiant:
                self._base_precision = self
            else:
                self._base_precision = ParentBigOh(self.base_ring())

    def __hash__(self):
        return id(self)

    def base_precision(self):
        return self._base_precision

    def precision(self):
        return self._precision

    def default_model(self):
        if self._default_mode is None:
            self.set_default_model()
        return self._default_model
    
    def set_default_model(self,model=None):
        if model is None:
            self._default_model = self._models.topological_sort()[-1]
        else:
            if self._models.has_vertex(model):
                self._default_model = model
            else:
                raise ValueError

    def add_model(self,model):
        from bigoh import BigOh
        if not isinstance(model,list):
            model = [model]
        for m in model:
            if not issubclass(m,BigOh):
                raise TypeError("A precision model must derive from BigOh but '%s' is not"%m)
            self._models.add_vertex(m)

    def delete_model(self,model):
        if isinstance(model,list):
            model = [model]
        for m in model:
            if self._models.has_vertex(m):
                self._models.delete_vertex(m)

    def update_model(self,old,new):
        from bigoh import BigOh
        if self._models.has_vertex(old):
            if not issubclass(new,BigOh):
                raise TypeError("A precision model must derive from BigOh but '%s' is not"%new)
            self._models.relabel({old:new})
        else:
            raise ValueError("Model '%m' does not exist"%old)

    def add_modelconversion(self,domain,codomain,constructor=None,safemode=False):
        if not self._models.has_vertex(domain):
            if safemode: return
            raise ValueError("Model '%s' does not exist"%domain)
        if not self._models.has_vertex(codomain):
            if safemode: return
            raise ValueError("Model '%s' does not exist"%codomain)
        path = self._models.shortest_path(codomain,domain)
        if len(path) > 0:
            raise ValueError("Adding this conversion creates a cycle")
        self._models.add_edge(domain,codomain,constructor)

    def delete_modelconversion(self,domain,codomain):
        if not self._models.has_vertex(domain):
            raise ValueError("Model '%s' does not exist"%domain)
        if not self._models.has_vertex(codomain):
            raise ValueError("Model '%s' does not exist"%codomain)
        if not self._models_has_edge(domain,codomain):
            raise ValueError("No conversion from %s to %s"%(domain,codomain))
        self._modelfs.delete_edge(domain,codomain)

    def uniformizer_name(self):
        return self._uniformizer

    def ambiant_space(self):
        return self._ambiant_space

    # ?!?
    def __call__(self,*args,**kwargs):
        return self._element_constructor_(*args,**kwargs)
    
    def _element_constructor_(self, *args, **kwargs):
        if kwargs.has_key('model'):
            del kwargs['model']
            return kwargs['model'](self, *args, **kwargs)
        if len(args) > 0:
            from precision.bigoh import BigOh, ExactBigOh
            arg = args[0]
            if isinstance(arg,BigOh) and arg.is_exact():
                return ExactBigOh(self)
            if self._models.has_vertex(arg.__class__) and arg.parent() is self:
                return arg
        if self._default_model is not None:
            try:
                return self._default_model(self,*args,**kwargs)
            except (AttributeError,ValueError,TypeError):
                pass
        models = self._models.topological_sort()
        models.reverse()
        for m in models:
            try:
                return m(self,*args,**kwargs)
            except (AttributeError,ValueError,TypeError,PrecisionError):
                pass
        raise PrecisionError("unable to create a BigOh object")

    def __repr__(self):
        return "Parent for BigOh in %s" % self._ambiant_space

    def models(self,graph=False):
        if graph:
            return self._models
        else:
            return self._models.vertices()

    def dimension(self):
        return self._ambiant_space.dimension()

    def indices_basis(self):
        return self._ambiant_space.indices_basis()
