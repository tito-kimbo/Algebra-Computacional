import random

from python_alcp.algorithms.discrete_log import discrete_log
from python_alcp.utils import externals

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
    

    elems = list(F.elements())
    elems.remove(F.zero)
    rho = random.choice(elems)
    while is_quadratic(rho):
        rho = random.choice(elems)

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

def modroot(r, e):
    """
        r-th root in a finite field.
        Adleman-Manders-Miller algorithm
        See https://arxiv.org/pdf/1111.4877.pdf
    """
    
    F = type(e)

    Z = externals.Z
    eea = externals.eea

    if e == F.zero:
        return [e]

    p = F.char()
    q = F.order()


    if (q-1) % r != 0:
        # Can't use the algorithm, fall back to discrete log
        
        return modroot_dl(r,e)

    else:
        if e**((q-1)//r) != F.one:
           res =  None 
        else:
            # Actual algorithm

            # Find an rth non-residue rho
            elems = list(F.elements())
            elems.remove(F.zero)
            rho = random.choice(elems)
            while (rho**((q-1)//r)).val.deg() == 0:# == F.one:
                rho = random.choice(elems)

            # Find s,t such that q-1 = s*r^t
            s = q-1
            t = 0
            while s % r == 0:
                s //= r
                t += 1

            if not eea(Z(r), Z(s))[0].is_unit():
                # Can't use algorithm, fall back to dl
                return modroot_dl(r,e)

            alpha = 0
            while (r*alpha - 1) % s != 0:
                alpha += 1

            #print(s,t)
            #print(alpha)

            #print(rho)

            a = rho**(s*r**(t-1))
            b = e**(r*alpha-1)
            c = rho**s
            h = 1
            for i in range(1,t):
                d = b**(r**(t-1-i))
                if d == F.one:
                    j = 0
                else:
                    #print(a,d)
                    j = -discrete_log(a,d)
                    #print(j)
                b = b*((c**r)**j)
                h = h*(c**j)
                c = c**r

            res = e**alpha * h

            if res is None:
                return None
            else:
                #print(rho)
                #print(r,s,t)
                a = rho**(s*(r**(t-1)))
                #print(a)
                #print([a**(i*r) for i in range(0,r)])
                return [res*(a**i) for i in range(0,r)]
    
def is_quadratic(e):
    """
        Test if an element of a finite field is a quadratic residue
    """
    return e == type(e).zero or e**((type(e).order()-1)//2) == type(e).one


def modroot_dl(r,e):
    F = type(e)

    if e == F.zero:
        return [e]

    p = F.char()
    q = F.order()

    Z = externals.Z
    eea = externals.eea
    x = F.generator()

    g, s, t = eea(Z(r), Z(q-1))

    if g.is_unit():
        res = e**(s / g).val

    else:
        dl = discrete_log(x,e)
        if dl is None or dl % g.val != 0:
            return None
        else:
            res = x**((dl//g.val) * s.val)

    return [res*(x**(i*(q-1)//g.val)) for i in range(g.val)]

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
