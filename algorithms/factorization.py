from collections import defaultdict
from algorithms.divisibility import gcd



def squarefree_decomposition(f):
    """
        f must be a monic polynomial in Fq[X]
    """

    def list_pow(d, p):
        # "Raises" a dict to a power as explained in Remark 2.3.8
        res = defaultdict(lambda: 1)
        for k,val in d.items():
            res[k*p] = val
        return res
        

    R = f.ring

    fp = f.der()

    p = R.char()

    if fp == R.zero:
        g = f.croot()
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

    g = gs[i].croot()

    L2 = list_pow(squarefree_decomposition(g), p)
    result = dict()
    for k,val in L:
        if k in L2:
            res[k] = val*L2[k]
        else:
            res[k] = val

    return result


def distinct_degree_decomposition(f):
    """
        f must be a squarefree polynomial in Fq[x]
    """
    
    R = f.ring
    CoefR = R.ring      # coefficient ring

    result = defaultdict(lambda: R.one)
    q = CoefR.order()
    d = 0
    g = [f]
    x = R.build([CoefR.zero, CoefR.one])
    h = [x]

    while d <= g[d].deg()/2 - 1:
        d += 1
        h.append(h[d-1]**q % g[d-1])
        fact = gcd(g[d-1], h[d] - x)
        if fact != R.one:
            result[d] = fact
        g.append(g[d-1] // fact)

    if g[d] != R.one:
        result[g[d].deg()] = g[d]

    return result
