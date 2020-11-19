from structures.rings import EuclideanDomain, IntegralDomain, UniqueFactorizationDomain
from functools import reduce


def eea(f: IntegralDomain.Element, g: IntegralDomain.Element):
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

def gcd(f: IntegralDomain.Element, g: IntegralDomain.Element, *args):
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

def gcd_dfu(*elems):
    """
        Greatest common divisor of a set of DFU elements
    """

    R = elems[0].ring
    if not isinstance(R, UniqueFactorizationDomain):
        raise ValueError("Arguments must be in UFD elements")

    def common_factors(a,b):
        common_keys = set(a.keys()).intersection(set(b.keys()))
        res = dict()
        for key in common_keys:
            res[key] = min(a[key], b[key])
        return res

    def multiply_all_factors(f):
        res = R.one
        for k,v in f.items():
            res *= k**v
        return res

    factors = reduce(common_factors, map(lambda x: x.factors(), elems))

    return multiply_all_factors(factors)
