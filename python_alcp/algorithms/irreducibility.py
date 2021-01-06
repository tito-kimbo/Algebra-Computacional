from python_alcp.utils import prime_factors, external, assuming
from python_alcp.algorithms.divisibility import gcd
from python_alcp.structures.ideals import GetQuotient, GetRingQuotient

@external
def rabin_test(pol):
    """
        Rabin's test of irreducibility for polynomials with coefficients in a finite field
    """

    PR = pol.ring       # polinomial ring
    R = PR.coefRing         # coefficient ring (field)
    assuming(R.is_finite())

    p = R.order()
    n = pol.deg()

    quotRing = GetRingQuotient(PR,(PR*pol))

    x = quotRing.build([R.zero,R.one])

    if (x**(p**n)-x) != quotRing.zero:
        return False

    factors = set(prime_factors(n).keys())
    for f in factors:
        ni = n//f
        if not gcd(pol, (x**(p**ni)-x).val).is_unit():
            return False

    return True
