import random
import itertools
from functools import reduce
from math import sqrt, ceil, log2


from python_alcp.examples.rings import Z
from python_alcp.algorithms.factorization import berlekamp_cantor_zassenhaus
from python_alcp.algorithms.divisibility import eea, gcd
from python_alcp.utils import assuming, primes


def hensel_lifting(f, g, h, s, t):
    """
        Input: polynomials in Z/(qZ)[x], q = p^r
        where f = gh, sg + th = 1
        Returns f', g', h', s' and t' in Z/(q^2)[x] lifted from g,h,s,t respectively
        which satisfy the same constraints
    """

    # Implements algorithm 2.4.6

    R = f.ring.coefRing
    
    q = R.ideal.generator

    newR = R.baseRing / (R.baseRing * (q**2))
    RX = newR[f.ring.var]


    # Move to the higher ring
    fp = RX(f)
    gp = RX(g)
    hp = RX(h)
    sp = RX(s)
    tp = RX(t)

    # Calculate lifted polynomials
    nabla = fp - gp * hp

    newg = gp * (RX.one + (sp*nabla // hp)) + tp*nabla
    newh = hp + (sp*nabla % hp)

    delta = sp*newg + tp*newh - RX.one
    news = sp - (sp * delta % newh)
    newt = (RX.one - delta) * tp - newg * (sp * delta // newh)

    # Check output meets requirements
    assuming(fp == newg * newh, "Hensel lifting failed")
    assuming(RX.one == news * newg + newt * newh, "Hensel lifting failed")

    return fp, newg, newh, news, newt




def recover_factors(f, gs):
    """
    """

    # Implements algorithm 2.4.10


    def represent(x, lb, ub):
        if x.val > ub:
            return x.val - x.ring.ideal.generator
        if x.val <= lb:
            return x.val + x.ring.ideal.generator
        return x.val

    R = gs[0].ring.coefRing
    RX = gs[0].ring
    ZX = Z[RX.var]

    order = R.order()
    lb = -(R.order() // 2)
    ub = (R.order() // 2)

    L = f.coefs[-1]
    h = ZX(f)
    d = 1
    I = list(range(len(gs)))

    result = []

    while 2*d <= len(I):
        P = list(itertools.combinations(I, d))
        while len(P) > 0 and 2*d <= len(I):
            S = random.choice(P)
            P.remove(S)
            g_ = reduce(RX.__mul__, [gs[i] for i in S], RX.build([L]))
            
            g = ZX.build([represent(v, lb, ub) for v in g_.coefs])

            if ZX.build([int(L)])*h % g == ZX.zero:
                gp = g.primitive_part()
                result.append(gp)
                h = h // gp
                for i in S:
                    I.remove(i)
                P = list(itertools.combinations(I, d))
                d += 1

    if h.deg() > 0:
        result.append(h)

    return result



def zx_factorization(f):
    """
    """
    # Implements Algorithm 2.4.3

    def norm(f):
        return sqrt(sum([int(a)**2 for a in f.coefs]))

    L = f.coefs[f.deg()]

    prim = primes()
    next(prim)
    p = next(prim)
    while True:
        RX = (Z/(p*Z))["x"]
        f_ = RX.build([a.val for a in f.coefs])
        if L.val % p != 0 and f_.der() != RX.zero and gcd(f_, f_.der()).is_unit():
            break
        p = next(prim)

    factors = berlekamp_cantor_zassenhaus(f_)

    nfac = len(factors)

    N = 1
    nor = norm(f)
    while p**N < (2**(f.deg()+1)) * nor:
        N += 1

    n_lifts = int(ceil(log2(N)))

    
    for i in range(n_lifts):
        # Hensel lifting from i to i+1

        RX = f_.ring
        RX.setRepr("reduced")
        newfactors = []
        nh = None
        nf = None

        for j in range(nfac-1):

            g = factors[j]
            h = reduce(RX.Element.__mul__, factors[j+1:], RX.one)

            gc, s, t = eea(g, h)

            if gc != gc.normal():
                s /= gc
                t /= gc

            assuming(gc.is_unit(), "Factorization failed")

            nf, ng, nh, ns, nt = hensel_lifting(f_, g, h, s, t)

            newfactors.append(ng)

        newfactors.append(nh)
        f_ = type(nf)(f_)
        factors = newfactors

    return recover_factors(f,factors)

