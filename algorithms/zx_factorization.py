from examples.rings import Z
from algorithms.factorization import berlekamp_cantor_zassenhaus
from algorithms.divisibility import eea, gcd
from utils import assuming
import random
import itertools
from functools import reduce
from math import sqrt

def _embed_polynomial(pol, target):
    coefs = target.ring
    res = target.build([coefs.build(a.val.val) for a in pol.coefs])
    return res

def hensel_lifting(f, g, h, s, t):
    """
        Input: polynomials in Z/(qZ)[x], q = p^r
        where f = gh, sg + th = 1
        Returns f', g', h', s' and t' in Z/(q^2)[x] lifted from g,h,s,t respectively
        which satisfy the same constraints
    """

    # Implements algorithm 2.4.6

    R = f.ring.ring
    
    q = R.ideal.generator

    newR = R.ring / (R.ring * (q**2))
    RX = newR["X"]


    # Move to the higher ring
    fp = _embed_polynomial(f, RX)
    gp = _embed_polynomial(g, RX)
    hp = _embed_polynomial(h, RX)
    sp = _embed_polynomial(s, RX)
    tp = _embed_polynomial(t, RX)

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

    R = f.ring.ring
    RX = f.ring
    ZX = Z["X"]

    L = f.coefs[f.deg()-1]
    h = _embed_polynomial(f, ZX)
    d = 1
    I = list(range(1,len(gs)+1))

    result = []

    while 2*d <= len(I):
        P = list(itertools.combinations(I, d))
        while len(P) > 0 and 2*d <= len(I):
            S = random.choice(P)
            P.remove(S)
            g_ = reduce(RX.Element.__mul__, [gs[i] for i in S], RX.build([L]))
            
            g = ZX.build([represent(v) for v in g_.coefs])

            if ZX.build([L.val.val])*h % g == ZX.zero:
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


curprimes = [2,3]
def primes():

    def _p():
        # Generator that returns primes
        
        count = 0
        while count < len(curprimes):
            yield curprimes[count]
            count += 1

        i = curprimes[-1] + 2
        while True:
            prime = True
            for p in curprimes:
                if i % p == 0:
                    prime = False
                    break

            if prime:
                curprimes.append(i)
                yield i

            i += 2

    return _p()
    

def zx_factorization(f):
    """
    """
    # Implements Algorithm 2.4.3

    def norm(f):
        return sqrt(sum([a.val**2 for a in f.coefs]))

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

    print(f_)
    print(RX)

    factors = berlekamp_cantor_zassenhaus(f_)

    nfac = len(factors)
    print(nfac)
    print(factors)

    N = 1
    nor = norm(f)
    while p**N < (2**(f.deg()+1)) * nor:
        N += 1

    
    for i in range(1,N):
        # Hensel lifting from i to i+1

        RX = f_.ring
        RX.setRepr("reduced")
        newfactors = []
        nh = None
        nf = None

        for j in range(nfac-1):

            g = factors[j]
            h = reduce(RX.Element.__mul__, factors[j+1:], RX.one)

            print(g, h)

            _, s, t = eea(g, h)

            print(s, t)
            print(_)

            assuming(_ == RX.one, "Factorization failed")

            nf, ng, nh, ns, nt = hensel_lifting(f_, g, h, s, t)

            newfactors.append(ng)

        newfactors.append(nh)
        f_ = _embed_polynomial(f_, nf.ring)
        factors = newfactors

    return recover_factors(factors)

