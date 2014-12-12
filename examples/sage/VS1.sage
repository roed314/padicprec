d = 5; dim = 3

K = Qp(2)
Kint = K.integer_ring()
V = K^d
Vint = Kint^d

E1 = VectorSpace(K,d)
E2 = VectorSpace(K,d)
gens1 = [ ]; gens2 = [ ]
for _ in range(dim):
    x = V.random_element()
    E1._add_in_place(x)
    gens1.append(x)
    x = V.random_element()
    E2._add_in_place(x)
    gens2.append(x)
W1 = V.subspace(gens1)
W2 = V.subspace(gens2)
