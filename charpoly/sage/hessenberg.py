class DynamicBasis:
    def __init__(self, base_ring, dimension):
        if base_ring not in DiscreteValuationRings():
            raise TypeError("%s is not a discrete valuation ring" % base_ring)
        self._base_ring = base_ring
        self._dimension = dimension
        self._space = base_ring**dimension
        self._vectors = [ ]

    def add_vector(self, v):
        v = self._space(v)
        coeffs = [ ]
        for (w, piv) in self._vectors:
            coeffs.append(v[piv])
            v -= v[piv]*w
        val = Infinity
        for i in range(self._dimension):
            n = v[i].valuation()
            if n < v[i].precision_absolute() and n < val:
                val = n
                piv = i
        if val is not Infinity:
            scalar = v[piv]
            v /= scalar
            self._vectors.append((v,piv))
            coeffs.append(scalar)
        else:
            v = None
        coeffs += [ self._base_ring(0) ] * (self._dimension - len(coeffs))
        return coeffs, v

    def basis(self):
        return [ v for (v,piv) in self._vectors ]


def hessenberg(M):
    if not M.is_square():
        raise ValueError
    dimension = M.nrows()
    ring = M.base_ring()
    space = ring**dimension
    
    basis = DynamicBasis(ring, dimension)
    hessenberg = [ ]
    i = 0; v = None
    while i < dimension:
        while v is None:
            v = space.random_element()
            _, v = basis.add_vector(v)
        while v is not None:
            v = M*v
            coeffs, v = basis.add_vector(v)
            hessenberg.append(coeffs)
            i += 1

    return matrix(hessenberg)
