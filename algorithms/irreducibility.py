def _factorize(n: int):
    # Crude factorization
    # TODO replace by sieve
    if n < 2:
        return []
    for i in range(2,n+1):
        if n % i == 0:
            return [i] + _factorize(n//i)

def rabin_test(pol):
    """
        Rabin's test of irreducibility for polynomials with coefficients in a finite field
    """

    from structures.ideals import RingQuotient

    PR = pol.ring       # polinomial ring
    R = PR.ring         # coefficient ring (field)
    assuming(R.is_finite())

    p = R.order()
    n = pol.deg()

    quotRing = RingQuotient(PR,(PR*pol))

    x = quotRing.build([R.zero,R.one])

    if (x**(p**n)-x) != quotRing.zero:
        return False

    factors = set(_factorize(n))
    for f in factors:
        ni = n//f
        if not gcd(pol, (x**(p**ni)-x).val).is_unit():
            return False

    return True
