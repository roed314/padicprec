from sage.rings.infinity import Infinity
from sage.matrix.constructor import matrix

def digits_of_precision(vectors):
    n = len(vectors[0])
    N = len(vectors)
    ring = vectors[0][0].parent()

    # Not diffused digits of precision
    not_diffused = n * [ Infinity ]
    for v in vectors:
        for i in range(n):
            val = v[i].valuation()
            if val < not_diffused[i]:
                not_diffused[i] = val
    not_diffused = sum(not_diffused)

    # All digits of precision
    all = 0
    M = matrix(vectors)
    for j in range(n):
        # Find pivot
        val = Infinity
        piv = 0
        for i in range(j,N):
            tmp = M[i,j].valuation()
            if tmp < val:
                val = tmp
                piv = i
        if M[piv,j] == 0:
            for i in range(j,N):
                print M[i,j], M[i,j].precision_absolute()
            raise RuntimeError

        all += val

        # Update matrix
        M.swap_rows(piv,j)
        for i in range(j+1,N):
            scalar = ring(- M[i,j] / M[j,j])
            M.add_multiple_of_row(i, j, scalar)

    diffused = all - not_diffused
    return all, diffused, not_diffused
