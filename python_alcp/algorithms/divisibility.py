from functools import reduce

from python_alcp.utils import external

@external
def eea(f, g):
    """ Extended Euclid Algorithm
        Given two elements of an euclidean domain, compute the gcd
        and the bezout identity elements.
        Returns (gcd, s, t) such that gcd = s*f + t*g
    """

    if f.ring != g.ring:
        raise ValueError("f and g must be elements of the same ring")

    R = f.ring

    r = [f, g]
    s = [R.one, R.zero]
    t = [R.zero, R.one]

    while r[-1] != R.zero:
        q = r[-2] // r[-1]
        r.append(r[-2] % r[-1])
        s.append(s[-2] - q*s[-1])
        t.append(t[-2] - q*t[-1])

    return (r[-2], s[-2], t[-2])

@external
def gcd(f, g, *args):
    """ Greatest common divisor, can be derived from EEA"""

    if f.ring != g.ring:
        raise ValueError("f and g must be elements of the same ring")
    R =  f.ring
    r = [f,g]
    while r[-1] != R.zero:
        r.append(r[-2] % r[-1])

    if len(args) > 0:
        return gcd(r[-2], *args)
    else:
        return r[-2]

#TODO: ensure both f and g are integers
#TODO: swap %, / and * to bitwise operations

@external
def binary_gcd(f, g, *args):
    """ Greatest common divisor, using the technique described in the 
        Pag 57, Exercise 4.1 of Shoup's book"""

    if f.ring != g.ring:
        raise ValueError("f and g must be elements of the same ring")
    R =  f.ring
    r = f
    r_ = g
    e = 0

    while r % 2 == 0 and r_ % 2 == 0:
        r = r / 2
        r_ = r_ / 2
        e = e + 1

    while True:
        while r % 2 == 0:
            r = r / 2
        while r_ % 2 == 0:
            r_ = r_ / 2

        if r_ < r:
            (r, r_) = (r_, r)

        r_ = r_ - r
        
        if r_ == R.zero:
            break

    d = 2**e*r

    if len(args) > 0:
        return binary_gcd(d, *args)
    else:
        return d
