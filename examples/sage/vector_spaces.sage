class VectorSpace:
    def __init__(self, base, dim):
        self._base = base
        self._ambiant_dim = dim
        self._gens = [ ]
        self._indices = range(dim)

    def __repr__(self):
        s = "Basis:\n"
        gens = MatrixSpace(self._base, self.dimension(), self.ambiant_dimension())(self._gens)
        s += gens.str()
        return s

    def copy(self):
        other = VectorSpace(self._base, self._ambiant_dim)
        other._gens = list(self._gens)
        other._indices = list(self._indices)
        return other

    def gens(self):
        return self._gens

    def ambiant_dimension(self):
        return self._ambiant_dim

    def dimension(self):
        return len(self._gens)

    def orthogonal(self):
        K = self._base
        dim = self.ambiant_dimension()
        d = dim - self.dimension()
        VS = FreeModule(K, dim)

        indices = list(self._indices)
        indices.reverse()

        gens = [ ]
        for i in range(d):
            index = indices[i]
            coords = [ K(0) ] * dim
            coords[index] = K(1)
            for j in range(d,dim):
                coords[indices[j]] = -self._gens[dim-1-j][index]
            vector = VS(coords)
            gens.append(vector)

        orth = VectorSpace(K, dim)
        orth._indices = indices
        orth._gens = gens
        return orth

    def __add__(self,other):
        ans = self.copy()
        ans._add_in_place(other.gens())
        return ans

    def _add_in_place(self, vectors):
        K = self._base
        gens = self._gens
        indices = self._indices
        dim = len(gens)
        ambiant_dim = self._ambiant_dim
        if not isinstance(vectors, list):
           vectors = [ vectors ]
        for vector in vectors:
            for i in range(dim):
                vector -= vector[indices[i]] * gens[i]
            val = Infinity; ii = 0
            for i in range(dim, ambiant_dim):
                coeff = vector[indices[i]]
                if coeff == 0: continue
                v = coeff.valuation()
                if val > v:
                    val = v; ii = i
            if val is Infinity: continue
            index = indices[ii]
            indices[dim], indices[ii] = index, indices[dim]
            vector /= vector[index]
            for i in range(dim):
                gens[i] -= gens[i][index] * vector
            gens.append(vector)
            dim += 1

    def intersection(self,other):
        V = self.orthogonal() + other.orthogonal()
        return V.orthogonal()

    def meet(self,other):
        return self.intersection(other)

    def apply(self,f):
        vectors = [ ]
        for vector in self._gens:
            vectors = f(vector)
        ans = VectorSpace(self._base, self._ambiant_dim)
        ans._add_in_place(vectors)
        return ans


d = 5
K = Qp(2)
Kint = K.integer_ring()
V = K^d
Vint = Kint^d

E1 = VectorSpace(K,d)
E2 = VectorSpace(K,d)
gens1 = [ ]; gens2 = [ ]
for _ in range(3):
    x = V.random_element()
    E1._add_in_place(x)
    gens1.append(x)
    x = V.random_element()
    E2._add_in_place(x)
    gens2.append(x)
W1 = V.subspace(gens1)
W2 = V.subspace(gens2)
