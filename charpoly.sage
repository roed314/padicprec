

def comatrix(A):
    M = A.parent()
    d = A.det()
    if not d:
        return A.parent().zero()
    n = A.ncols()
    R = M.base_ring()
    K = R.fraction_field()
    B = A.change_ring(K).augment(M.one() * d)
    return B.echelon_form().change_ring(R).matrix_from_columns(range(n,2*n))

def prec_gens(A):
    """
    INPUT:

    - ``A`` -- a square matrix over Z

    OUTPUT:

    - the index of the precision increase inside Z^n
    - a sorted list of generators for the precision increase
    """
    R = A.parent().base_ring()
    Rx = PolynomialRing(R,'x')
    x = Rx.gen()
    n = A.ncols()
    B = x*identity_matrix(n) - A
    zero = R(0)
    C = comatrix(B)
    def padded(L):
        if len(L) < n:
            return L + [zero] * (n - len(L))
        else:
            return L[:n]
    gen_mat = matrix(R,[padded(c.list()) for c in C.list()]).transpose().change_ring(R)
    D, U, V = gen_mat.smith_form()
    rgen_mat = (~U)*D
    index = prod([D[i,i] for i in range(n)])
    return index, sorted([Rx(list(rgen_mat.column(c))) for c in range(n)])

from collections import defaultdict

class NonmaximalExamples(SageObject):
    def __init__(self):
        self.trials = defaultdict(int)
        self.examples = defaultdict(lambda: defaultdict(list))

    def add_examples(self, dim, N=1024, num_examples=10, p=None, index_sought=None):
        """
        INPUT:

        - ``dim`` -- the dimension of the matrices to search for
        - ``N`` -- bound up to which integers are sampled uniformly in creating random matrices
        - ``num_examples`` -- the number of new examples to find
        - ``p`` -- a prime, given when searching for a particular index valuation
        """
        trials = 0r
        found = 0r
        examples = self.examples[dim,N]
        while found < num_examples:
            A = random_matrix(ZZ,dim,x=N)
            if A.det() == 0: continue
            index, G = prec_gens(A)
            if index != 1r:
                examples[index].append((A,G))
                found += 1r
                if index_sought is not None:
                    if index.valuation(p) >= index_sought.valuation(p):
                        return (A, G)
            trials += 1r
        self.trials[dim,N] += trials
    def add_examples_manydim(self, dims, iterations=500):
        for t in range(iterations):
            for d in dims:
                self.add_examples(d)

    def add_smith_examples(self, dim, p, max_val=40, num_examples=10):
        trials = 0r
        found = 0r
        examples = self.examples[dim,0]
        while found < num_examples:
            vals = []
            while len(vals) < dim:
                new_val = ZZ.random_element()
                if 0 <= new_val <= max_val:
                    vals.append(p^new_val)
            D = diagonal_matrix(sorted(vals))
            U = random_matrix(ZZ,dim,algorithm="unimodular")
            V = random_matrix(ZZ,dim,algorithm="unimodular")
            A = U * D * V
            index, G = prec_gens(A)
            if index != 1r:
                examples[index].append((A,G))
                found += 1r
            trials += 1r
        self.trials[dim,0] += trials

    @staticmethod
    def _match_dim(dim, d):
        return dim is None or d == dim
    @staticmethod
    def _match_N(p, N, strict):
        return N != 0 and (not strict or N.is_power_of(p))

    def index_stats(self, p, dim=None, strict=False, percentage=False):
        def match(d, N):
            return self._match_dim(dim, d) and self._match_N(p, N, strict)
        index_vals = defaultdict(int)
        total = 0r
        for key, examples in self.examples.iteritems():
            if match(*key):
                T = self.trials[key]
                total += T
                for index in examples.iterkeys():
                    index_vals[index.valuation(p)] += 1r
                    T -= 1
                # The unstored examples have index 1
                index_vals[0] += T
        if not index_vals: return []
        max_val = max(index_vals.keys())
        if percentage:
            total = RR(total)
            return [RR(index_vals[v]) / total for v in range(max_val+1)]
        else:
            return [index_vals[v] for v in range(max_val+1)]

    @staticmethod
    def _match_diag(need_diag, example):
        if need_diag is None: return True
        actually_diag = all(g.is_term() for g in example[1])
        return (need_diag and actually_diag) or (not need_diag and not actually_diag)
    @staticmethod
    def _match_index(p, index_required, actual_index):
        if index_required is None: return True
        return index_required.valuation(p) == actual_index.valuation(p)

    def find_example(self, p, dim, index=None, diagonal=None, all=False, k=10, max_extra_examples=1000):
        L = []
        for key, examples in self.examples.iteritems():
            if self._match_dim(dim, key[0]):
                for ind, exs in examples.iteritems():
                    if self._match_index(p, index, ind):
                        for ex in exs:
                            if self._match_diag(diagonal,ex):
                                if all:
                                    L.append(ex)
                                else:
                                    return ex
        if not all:
            if dim is not None:
                ex = self.add_examples(dim, p^k, max_extra_examples, p, index)
                if ex is not None:
                    return ex
            raise RuntimeError("No example found")
        return L

    def all_diagonal_examples(self):
        L = []
        for key, examples in self.examples.iteritems():
            for ind, exs in examples.iteritems():
                for ex in exs:
                    if self._match_diag(True, ex):
                        L.append(ex)
        return L

    def all_nondiagonal_examples(self):
        L = []
        for key, examples in self.examples.iteritems():
            for ind, exs in examples.iteritems():
                for ex in exs:
                    if self._match_diag(False, ex):
                        L.append(ex)
        return L

    def all_smith_examples(self):
        L = []
        for key, examples in self.examples.iteritems():
            if key[1] == 0:
                for ind, exs in examples.iteritems():
                    L.extend(exs)
        return L

def update(NME):
    NME2 = NonmaximalExamples()
    NME2.examples = NME.examples
    NME2.trials = NME.trials
    return NME2

def minval(A, p):
    return min([c.valuation(p) for c in A.list()])

def convexify(L):
    P = Polyhedron(vertices = list(enumerate(L)), rays=[(0,1)])
    verts = sorted(P.vertices_list())
    heights = [QQ(verts[0][1])]
    for i in range(len(verts)-1):
        cur = verts[i]
        next = verts[i+1]
        slope = (next[1] - cur[1])/(next[0] - cur[0])
        for x in range(cur[0],next[0]):
            heights.append(heights[-1] + slope)
    return heights

def NP(A, p):
    f = A.charpoly()
    return convexify([a.valuation(p) for a in f])

def HP(A, p):
    D, U, V = A.smith_form()
    n = min(D.nrows(), D.ncols())
    vals = sorted([D[i,i].valuation(p) for i in range(n)])
    S = sum(vals)
    HP = [S]
    for v in reversed(vals):
        HP.append(HP[-1] - v)
    return HP

def CP(G, p):
    vals = [min([g[i].valuation(p) for g in G if g[i]]) for i in range(len(G))]
    return convexify(vals)

def test_between(L, p):
    for A, G in L:
        np = NP(A, p)
        hp = HP(A, p)
        cp = CP(G, p)
        for i in range(len(cp)):
            if cp[i] < hp[i+1] or cp[i] > np[i+1]:
                print A
                print G
                raise RuntimeError

def HPexceed(L, p):
    exceed = defaultdict(int)
    for A, G in L:
        hp = HP(A, p)
        cp = CP(G, p)
        diff = tuple(cp[i] - hp[i+1] for i in range(len(cp)))
        exceed[diff] += 1
    return exceed

def HPexceed_display(L, p):
    exceed = HPexceed(L, p)
    excesses = sorted(exceed.keys(), cmp=lenlex)
    for key in excesses:
        print key, exceed[key]

def lenlex(x, y):
    # For sorting difference tuples
    c = cmp(len(x),len(y))
    if c: return c
    return cmp(x, y)

def print_polygons(ex, p):
    A, G = ex
    print "NP: ", NP(A,p)
    print "HP: ", HP(A,p)
    print "CP:    ", CP(G,p)
