from collections import defaultdict
from algorithms.divisibility import gcd
import random
from math import log2



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


def distinct_degree_factorization(f):
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


def berlekamp_splitting(f, hs):
    """
        f must be squarefree and monic
        hs must be a basis of ker(phi_f))
    """
    # Implements algorithm 2.3.21

    # R = ring of coefficients, RX = polynomial ring
    R = f.ring.ring
    RX = f.ring

    # s is the dim of Ker(phi) and must be
    # the number of factors, per Lemma 2.3.19
    s = len(hs)
    if s == 1:
        return f

    
    # Compute all elements of R = Fq.
    alpha = R.generator()
    R_elems = [R.zero] + [alpha**i for i in range(R.order()-1)]

    # Loop through the basis and R_elems until we find a nontrivial factor
    i = 1
    while i < s:
        for a in R_elems:
            g = gcd(f, hs[i]-RX.build([a]))
            if g != RX.one and g != f:
                return g
        i += 1

    raise ValueError("You broke mathematics")

def ker_phi_basis(f):
    """
        Computes a basis for Ker(Phi_f),
        where is the endomorphism (Fr - id)
        in Fq[x] / <f>
    """

    # Implements Remark 2.3.22
    # Computes the matrix of Phi_f and then applies
    # gaussian elimination to find independent zeroes of Phi_f


    # R = ring of coefficients, RX = polynomial ring
    R = f.ring.ring
    RX = f.ring

    d = f.deg()
    q = R.order()

    x = RX.build([R.zero,R.one])

    # B is a basis of R, M is the matrix of Phi_f in that basis
    # Initialize B as the standard basis and M accordingly
    B = [[R.zero]*i + [R.one] + [R.zero]*(d-1-i) for i in range(d)]
    M = [((x**(i*q) - x**i) % f).coefs for i in range(d)]

    # M is build as a vector of polynomials, so we need
    # to add leading zero coefficients if the degree is
    # not high enough
    for i in range(d):
        if len(M[i]) < d:
            M[i] += [R.zero] * (d-len(M[i]))

    # Modifies B and M in-place
    _gaussian_elimination(B,M,R)

    # If M[i] = 0, then Phi_f(B[i]) = 0
    zeros = [i for i in range(d) if M[i] == [R.zero]*d]
    return [RX.build(B[i]) for i in zeros]


def berlekamp_factorization(f):
    """
        Berlekamp factorization algorithm (BFA)
    """
    # Implements algorithm 2.3.24

    # Compute a basis of Ker(Phi)
    basis = ker_phi_basis(f)

    # s must be the number of factors
    s = len(basis)

    # We will store found factors in factors,
    # and irreducible ones in irred
    # TODO Better to use sets, but hash needed
    factors = [f]
    irred = []
    
    # While we have less than s factors, attemp to reduce one in factors
    while len(factors) + len(irred) < s:
        g = random.choice(factors)
        h = berlekamp_splitting(g, basis)
        factors.remove(g)
        if g != h:
            factors.append(h)
            factors.append(g // h)
        else:
            # If the splitting found no nontrivial factors, g must be irreducible
            irred.append(g)   

    return factors + irred


def berlekamp_cantor_zassenhaus(f):
    """
        Berlekamp / Cantor / Zassenhaus factorization algorithm (BCZ)
    """

    # Implements Algorithm 2.3.25

    # R = ring of coefficients, RX = polynomial ring
    R = f.ring.ring
    RX = f.ring

    q = R.order()
    
    # There is one operation which is different in characteristic 2
    if q % 2 == 0:
        r = int(log2(q))
        def op(x):
            return sum([x**(2*i) for i in range(r-1)], x.ring.zero)
    else:
        def op(x):
            return x**((q-1)//2) - x.ring.one


    # Compute a basis of Ker(Phi_f)
    hs = ker_phi_basis(f)

    # s must be the number of factors
    s = len(hs)

    # Store factors in result
    result = [f]

    alpha = R.generator()

    # Look for factors until we have found all
    while len(result) < s:
        g = random.choice(result)
        while g.deg() <= 1:
            g = random.choice(result)

        # Construct a random element of Ker(Phi)
        cs = [alpha**random.randint(0,q) for i in range(s)]
        h = sum([RX.build([a])*b for a,b in zip(cs,hs)], RX.zero)

        w = gcd(g, op(h))

        if w != RX.one and w != g:
            # w is a nontrivial factor of g
            result.remove(g)
            result.append(w)
            result.append(g//w)

    return result




def _gaussian_elimination(B, M, R):
    """
        Note: modifies B and M in place
    """

    l = len(M)

    for i in range(l):
        # find pivot
        pivots = [j for j in range(i,l) if M[j][i] != R.zero]
        if len(pivots) == 0:
            continue
        else:
            p = pivots[0]
            # swap rows
            M[p], M[i] = M[i], M[p]
            B[p], B[i] = B[i], B[p]

            # eliminate
            for j in pivots[1:]:
                c = M[j][i] // M[i][i]
                for k in range(i):
                    B[j][k] -= c*B[i][k]
                for k in range(i, l):
                    M[j][k] -= c*M[i][k]
                    B[j][k] -= c*B[i][k]
