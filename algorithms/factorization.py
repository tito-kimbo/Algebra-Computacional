from collections import defaultdict
from algorithms.divisibility import gcd



def squarefree_decomposition(f):
    """
        f must be a monic polynomial in Fq[X]
    """

    def list_pow(d, p):
        # "Raises" a dict to a power as explained in Remark 2.3.8
        res = defaultdict(lambda: 1)
        for k,val in d:
            res[k*p] = val
        return res
        

    R = f.ring

    fp = f.der()

    p = R.char()

    if fp == R.zero:
        g = fp.root(p)
        return list_pow(squarefree_decomposition(g), p)

    i = 0
    L = dict()
    gs = [gcd(f,fp)]
    ws = [f // gs[0]]

    while ws[i] != R.one:
        i = i+1
        ws.append(gcd(gs[i-1],ws[i-1]))
        gs.append(gs[i-1] // ws[i])
        L[i] = ws[i-1] // ws[i]

    if gs[i] == R.one:
        return L

    g = gs[i].root(p)

    L2 = list_pow(squarefree_decomposition(g), p)
    result = dict()
    for k,val in L:
        if k in L2:
            res[k] = val*L2[k]
        else:
            res[k] = val

    return result
