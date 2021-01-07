def modsqrt(e):
    """
        Square root in a finite field.
        Adleman-Manders-Miller algorithm
        See https://arxiv.org/pdf/1111.4877.pdf
    """
    
    if not is_quadratic(e):
        return None

    F = type(e)

    if e == F.zero:
        return [e]

    p = F.char()
    q = F.order()
    
    rho = None
    for elem in F.elements():
        if elem != F.zero and not is_quadratic(elem):
            rho = elem
            break

    s = q-1
    t = 0
    while s % 2 == 0:
        s //= 2
        t += 1

    a = rho**s
    b = e**s
    h = 1
    for i in range(1,t):
        d = b**(2**(t-1-i))
        k = 0 if d == F.one else 1
        b = b*((a**2)**k)
        h = h*(a**k)
        a = a**2

    res = e**((s+1)//2)*h

    if res**2 != e:
        raise AssertionError("Could not find sqrt of {e} in {type(e)}. Is the field correctly constructed?")
    return [res,-res]


def is_quadratic(e):
    """
        Test if an element of a finite field is a quadratic residue
    """
    return e == type(e).zero or e**((type(e).order()-1)//2) == type(e).one


def shifting_root(n, e, b):
    """
        Finds the nth root of e in base b using the shifting algorithm
        e must be represented as a list of digits in base b
    """

    def sum_base(e,b):
        return sum([b**i * x for i,x in enumerate(e[::-1])], b.ring.zero)

    R = b.ring

    x = R.zero
    y = R.zero
    r = R.zero

    mod = len(e) % n
    if mod > 0:
        alpha = sum_base(e[0:mod], b)
        i = mod
    else:
        alpha = sum_base(e[0:n], b)
        i = n

    while i <= len(e):

        x = b**n * x + alpha
        
        beta = R.zero
        while (b*y + beta)**n < x:
            beta = beta + R.one

        y = b*y + beta
        r = x - y**n

        i += n
        alpha = sum_base(e[i-n:i], b)

    return y
