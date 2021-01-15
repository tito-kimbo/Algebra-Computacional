from python_alcp.structures.multipoly import *

"""
Returns a tuple (is_interreduced,i,j). The first item indicates whether the given set of polynomials is interreduced. 
If it isn't, the polynomials at positions i and j verify that LM(F[i]) divides 0LM(F[j]).

:F list of polynomials over the same polynomial ring.
"""
def is_interreduced(F):
    # Does not contain 0
    for f in F:
        if f == f.ring.zero:
            return False
    for i in range(len(F)):
        for j in range(i+1,len(F)):
            if divides_monomial(F[i].lm(),F[j].lm()):
                return False,i,j    
            elif divides_monomial(F[j].lm(),F[i].lm()):
                return False,j,i
    return True,-1,-1

"""
Returns an interreduced list of polynomials w.r.t the given one.

:F list of polynomials over the same polynomial ring
"""
def interreduce(F):
    A = set(F)
    if len(A)==0:
        return list(A)
    
    # remove 0
    zero = type(next(iter(A)))(0)
    if zero in A:
        A.remove(zero)
    lA = list(A)
    stop,i,j = is_interreduced(lA)
    while not stop:
        # choose f != g in A with lm(f)|lm(g) -> lm(A[i]) | lm(A[j])
        # update g
        f,g = lA[i],lA[j]
        A.remove(g)
        g = g-term_to_poly(tuple_div(g.lt(),f.lt()),type(zero))*f
        if g != zero:
            A.add(g)
        lA = list(A)
        stop,i,j = is_interreduced(lA)
    return lA

"""
Returns the LCM of 2 monomials.

:m1 monomial
:m2 monomial
"""
def monomial_lcm(m1,m2):
    deg = tuple([max(m1.deg[i],m2.deg[i]) for i in range(len(m1.deg))])
    return Monomial(deg,m1.vars)

"""
Returns the s-polynomial with respect to 2 given polynomials.

:f  polynomial
:g  polynomial
"""    
def s_poly(f,g):
    x_gamma = monomial_lcm(f.lm(),g.lm())
    m1 = Monomial(deg_diff(x_gamma,f.lm()), f.lm().vars)
    m2 = Monomial(deg_diff(x_gamma,g.lm()), f.lm().vars)
    t1 = term_to_poly((m1,g.lc()),type(f))
    t2 = term_to_poly((m2,f.lc()),type(f))
    return t1*f-t2*g

"""
Returns a Gröbner basis for the ideal generated by F.

:F list of polynomials over the same polynomial ring
"""
def buchberger(F):
    G = set(F)
    P = set([(F[i],F[j]) for i in range(len(F)) for j in range(len(F)) if j>i])
    while len(P)>0:
        p = next(iter(P)) # choose an element of P
        f,g = p
        # update P
        P.remove(p)
        # Set h
        h = div_poly_RNF(s_poly(f,g),list(G))[1]
        if h != h.ring.zero:
            # update P and G
            for f in G:
                P.add((h,f))
            G.add(h)
    return list(G)