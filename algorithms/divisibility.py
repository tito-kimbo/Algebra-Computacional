from structures.rings import EuclideanDomain, IntegralDomain


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
        q = r[-2] / r[-1]
        r.append(r[-2] % r[-1])
        s.append(s[-2] - q*s[-1])
        t.append(t[-2] - q*t[-1])

    return (r[-2], s[-2], t[-2])

def gcd(f: IntegralDomain.Element, g: IntegralDomain.Element):
    """ Greatest common divisor, can be derived from EEA"""
    if f.ring != g.ring:
        raise ValueError("f and g must be elements of the same ring")
    R =  f.ring
    r = [f,g]
    while r[-1] != R.zero:
        q = r[-2] / r[-1]
        r.append(r[-2] % r[-1])
    return r[-2]

def eea_test():

    from examples.rings import Z
    from algorithms.eea import eea

    z = Z()


    f = z.SymbolicInteger(68)
    g = z.SymbolicInteger(129)
    gcd,s,t = eea(f, g)
    assert(gcd.val == 1)
    assert(gcd == s*g + t*f)
    print(f"{str(gcd)} is the gcd of {str(f)} and {str(g)}")
    print(f"{str(gcd)} = {str(t)}*{str(f)} + {str(s)}*{str(g)}")

    f = z.SymbolicInteger(128)
    g = z.SymbolicInteger(24)
    gcd,s,t = eea(f, g)
    assert(gcd.val == 8)
    assert(gcd == s*g + t*f)
    print(f"{str(gcd)} is the gcd of {str(f)} and {str(g)}")
    print(f"{str(gcd)} = {str(t)}*{str(f)} + {str(s)}*{str(g)}")
