from itertools import zip_longest

from python_alcp.structures.rings import (
        Ring,
        IntegralDomain,
        EuclideanDomain,
        RingElement,
        EuclideanDomainElement
)
from python_alcp.utils import (
        assuming,
        print_superscript,
        external,
        externals,
        prime_factors
)

VARS = ["X","Y","Z","T","U","V"] # Could be extended arbitrarily with sub indexing

def addAll(l,neutral):
    res = neutral
    for x in l:
        res = res+x
    return res


# Polynomials over integral domains
class PolynomialRing(Ring):

    def __eq__(cls, other):
        return hasattr(other, "coefRing") and cls.coefRing == other.coefRing

    def char(cls):
        return cls.coefRing.char()

    def order(cls):
        return 0

    def setRepr(cls, rep):
        cls.repr = rep
        if rep == "reduced":
            cls.coefRing.setRepr(rep)

    def generators(cls):
        return {cls(0,1)}.union({cls(c) for c in cls.coefRing.generators()})

    def units(cls):
        return {cls(u) for u in cls.coefRing.units()}

    def is_polynomial(cls):
        return True


class PolynomialED(PolynomialRing, EuclideanDomain):
    """
        A polynomial ring which is an euclidean domain, or equivalently, whose ring of coefficients is a field
    """

    def phi(cls, element):
        return element.deg()

    def numphi(cls, n: int):
        return cls.char() ** n




class PolynomialRingElement(RingElement):

    def __init__(self,*val):
        R = type(self).coefRing

        if len(val) == 1:
            if hasattr(val[0], "__iter__"):
                val = val[0]
            elif hasattr(val[0], "coefs"):
                val = val[0].coefs
            while hasattr(val, "val") and hasattr(val.val, "__iter__"):
                val = val.val

        cs = [c if isinstance(c,R) else R(c) for c in val]

        while len(cs) > 0 and cs[-1] == R.zero:
            cs.pop()

        self.val = tuple(cs)
        self.coefs = tuple(cs)
    
    def deg(self):
        return len(self.val)-1
    
    def __add__(self,other):
        if not hasattr(other, "val") or not hasattr(other.val, "__iter__"):
            other = type(self)(other)
        return type(self)([x+y for x,y in zip_longest(self.val,other.val,fillvalue=type(self).coefRing.zero)])

    def __sub__(self,other):
        if not hasattr(other, "val") or not hasattr(other.val, "__iter__"):
            other = type(self)(other)
        return type(self)([x-y for x,y in zip_longest(self.val,other.val,fillvalue=type(self).coefRing.zero)])
    
    # Convolutional product
    def inner_mul(self,other):
        if not hasattr(other, "val") or not hasattr(other.val, "__iter__"):
            other = type(self)(other)
        R = type(self).coefRing
        newcs = [R.zero] * (self.deg() + other.deg() + 1)
        for i,v in enumerate(self.val):
            for j,w in enumerate(other.val):
                newcs[i+j] += v*w

        return type(self)(newcs)

    def __neg__(self):
        return type(self)([-x for x in self.val])

    def __floordiv__(self, other):
        if not hasattr(other, "val") or not hasattr(other.val, "__iter__"):
            other = type(self)(other)
        quot, rem = polynomial_division(self, other)
        return quot

    def __mod__(self, other):
        if not hasattr(other, "val") or not hasattr(other.val, "__iter__"):
            other = type(self)(other)
        quot, rem = polynomial_division(self, other)
        return rem

                
    def __eq__(self,other):
        return hasattr(other, "val") and self.val == other.val
    
    def __str__(self):
        if type(self).repr == "reduced":

            if self == type(self).zero:
                return "0"

            reps = [f"{'' if a == a.ring.one and i > 0 else repr(a)}{type(self).var if i > 0 else ''}{print_superscript(i) if i > 1 else ''}" for i,a in enumerate(self.val) if a != a.ring.zero]

            return " + ".join(reversed(reps))

        else :# type(self).repr == None:

            if self == type(self).zero:
                return "(0)"

            reps = [f"{repr(a)}{'*' if i > 0 else ''}{type(self).var if i > 0 else ''}{'^' if i > 1 else ''}{i if i > 1 else ''}" for i,a in enumerate(self.val) if a != a.ring.zero]

            return "(" + " + ".join(reversed(reps)) + ")"


    def __lt__(self, other):
        if self.deg() >= other.deg():
            return False
        else:
            for i in range(1,self.deg()+1):
                if self.val[-i] >= self.val[-i]:
                    return False
            return True

    def inverse(self):
        assuming(self.is_unit(), "Can't invert non-unit")
        return type(self)(self.val[0].inverse())

    def der(self):
        # derivative
        new_val = []
        for i,c in enumerate(self.val[1:]):
            if c != type(self).zero:
                new_val.append((i+1)*c)
        return type(self)(new_val)


    def is_unit(self):
        return self.deg() == 0 and self.val[0].is_unit()

    def normal(self):
        # First we take the normal form of the leading coefficient
        nf = self.val[-1].normal()

        # If both are the same, the polynomial is already in it's normal form
        if nf == self.val[-1]:
            return self
        
        # If not, divide by leading coef
        else:
            div = (self.val[-1] / nf)
            poly = type(self)([c / div for c in self.val])
            return poly


    def content(self):
        if len(self.val) == 0:
            return type(self).coefRing.zero
        elif len(self.val) == 1:
            return self.val[0]
        else:
            return externals.gcd(*self.val)

    def primitive_part(self):
        cont = self.content()
        cf = list(map(lambda x: x // cont, self.val))
        return type(self)(cf)

    def is_prime(self):
        if type(self).coefRing.is_finite():
            return externals.rabin_test(self)
        else:
            raise NotImplementedError()

    def __hash__(self):
        return hash((type(self).__name__, self.val))



class PolynomialEDElement(PolynomialRingElement, EuclideanDomainElement):

    def __floordiv__(self, other):
        if not hasattr(other, "val") or not hasattr(other.val, "__iter__"):
            other = type(self)(other)
        quot, rem = polynomial_division(self, other)
        return quot

    def __mod__(self, other):
        if not hasattr(other, "val") or not hasattr(other.val, "__iter__"):
            other = type(self)(other)
        quot, rem = polynomial_division(self, other)
        return rem

    def is_unit(self):
        return self.deg() == 0

    def factors(self):
        raise NotImplementedError()

    def croot(self):
        """ cth root of the element, where c is the characteristic of the ring """
        R = type(self)
        p = R.char()
        newc = []
        k = 0

        for i,c in enumerate(self.val):
            if i == k*p:
                newc.append(c.croot())
                k += 1

        return R.build(newc)

    def is_primitive_field(self):
        if not self.is_prime():
            return False
        if not type(self).coefRing.is_finite():
            return False

        F = externals.GetFieldQuotient(type(self), type(self)*self)
        O = F.order()-1
        factors = prime_factors(O)
        x = F.generator()

        for f in factors.keys():
            if x**(O//f) == F.one:
                return False

        return True



def polynomial_division(a, b, pseudo = False):
    """ 
        Computes a // b, where a and b are polynomials.
        Returns (quotient, remainder)
        Long division algorithm
    """
    
    R = type(a).coefRing

    if b == type(b).zero:
        raise ValueError("Can't divide polynomial by 0")

    if a.deg() < b.deg():
        return (type(a).zero, a)

    delta = a.deg() - b.deg() + 1
    rem = list(a.val)
    div = list(b.val)
    
    quot = [R.zero] * (delta)

    lcdiv = div[-1]
    while len(rem) >= len(div):

        lc = rem[-1]
        degdif = len(rem) - len(div)
        nquot = R(lc / lcdiv) if lcdiv.is_unit() else R(lc // lcdiv)

        if lc != lcdiv*nquot:
            # Division is undefined, try pseudodivision
            if not pseudo:
                raise ValueError(f"Division undefined for {a} and {b} in {R}")
            return polynomial_division(a*(b.val[-1]**delta),b, pseudo)

        for i,v in enumerate(div):
            rem[i+degdif] -= v*nquot

        while len(rem) > 0 and rem[-1] == R.zero:
            rem.pop()
        quot[degdif] += nquot

    return type(a)(quot), type(a)(rem)

@external
def GetPolynomials(ring, var = None):
    if var is None:
        if hasattr(ring, "chain"):
            var = VARS[ring.chain+1]
            chain = ring.chain+1
        else:
            var = VARS[0]
            chain = 0
    else:
        chain = 0
    attrs = {"coefRing": ring, "var": var, "chain": chain}
    if ring.is_field():
        return PolynomialED(f"{ring}[{var}]", (PolynomialEDElement,), attrs)
    else:
        return PolynomialRing(f"{ring}[{var}]", (PolynomialRingElement,), attrs)
