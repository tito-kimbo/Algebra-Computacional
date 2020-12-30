from structures.rings import Ring, IntegralDomain, EuclideanDomain
from structures.fields import Field
from itertools import zip_longest
from algorithms.divisibility import gcd
from utils import assuming, print_superscript

VARS = ["X","Y","Z","T","U","V"] # Could be extended arbitrarily with sub indexing

def addAll(l,neutral):
    res = neutral
    for x in l:
        res = res+x
    return res


# Polynomials over integral domains, WIP
class PolynomialRing(Ring):

    # Univariate polynomials over a Ring - dynamically linked to the current ring
    class Element(Ring.Element):

        def __init__(self,coefs):
            R = self.ring.ring
            for c in coefs:
                assuming(c.ring == R, "coefficients must be ring elements")

            while len(coefs) > 0 and coefs[-1] == R.zero:
                coefs.pop()

            self.coefs = coefs
        
        def deg(self):
            return len(self.coefs)-1
        
        def __add__(self,other):
            return self.__class__([x+y for x,y in zip_longest(self.coefs,other.coefs,fillvalue=self.ring.ring.zero)])

        def __sub__(self,other):
            return self.__class__([x-y for x,y in zip_longest(self.coefs,other.coefs,fillvalue=self.ring.ring.zero)])
        
        # Convolutional product
        def inner_mul(self,other):
            D = self.deg()+other.deg()+1
            c,o = self.coefs[:],other.coefs[:] # shallow copies to avoid issues with original
            ring = self.ring.ring
            while len(c) < D:
                c.append(self.ring.ring.zero)
            coefs = list(zip_longest(c,o,fillvalue=ring.zero))
            aux = []
            for i in range(D):
                aux.append(addAll([coefs[k][0]*coefs[i-k][1] for k in range(i+1)],ring.zero))
            while len(aux)!=0 and aux[-1] == ring.zero:
                aux.pop(-1)
            return self.__class__(aux)

        def __neg__(self):
            return self.__class__(list(map(lambda x: -x, self.coefs)))

                    
        def __eq__(self,other):
            return isinstance(other, Ring.Element) and self.ring == other.ring and self.coefs == other.coefs
        
        def __str__(self):
            if self.ring.repr == "reduced":
                if len(self.coefs) == 0:
                    return "0"

                total_s = []

                if self.coefs[0] != self.ring.ring.zero:
                    s = str(self.coefs[0])
                    if " " in s:
                        s = "("+s+")"
                    total_s.append(s)


                for i in range(1,self.deg()+1):
                    s = ""
                    if self.coefs[i] != self.ring.ring.zero:
                        c = str(self.coefs[i])
                        if " " in c:
                            c = "("+c+")"
                        if self.coefs[i] == self.ring.ring.one:
                            c = ""
                        s += c  + self.ring.var + print_superscript(i)
                        total_s.append(s)

                return " + ".join(total_s)

            else :# self.ring.repr == None:
                if len(self.coefs) == 0:
                    return "(0)"
                s = "(" + str(self.coefs[0])
                for i in range(1,self.deg()+1):
                    if self.coefs[i] != self.ring.ring.zero:
                        c = str(self.coefs[i]) 
                        if self.coefs[i] == self.ring.ring.one:
                            c = ""
                        s += " + " + c + "*" + self.ring.var + "^" + str(i)
                s += ")"
                return s


        def der(self):
            # derivative
            new_coefs = []
            for i,c in enumerate(self.coefs[1:]):
                if c != self.ring.zero:
                    new_coefs.append((i+1)*c)
            return self.__class__(new_coefs)


        def is_unit(self):
            raise NotImplementedError()

        def normal(self):
            #First we take the normal form of the leading coefficient
            nf = self.coefs[-1].normal()

            #If both are the same, the polynomial is already in it's normal form
            if nf == self.coefs[-1]:
                return self
            
            else:
                poly = self * self.ring.build([nf / self.coefs[-1]]) 
                return poly

        def content(self):
            return gcd(*self.coefs)

        def primitive_part(self):
            cont = self.content()
            cf = list(map(lambda x: x / self.content(), self.coefs))
            return self.ring.build(cf)

    def __str__(self):
        return f"{self.ring}[{self.var}]"

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return self.ring == other.ring

    def char(self):
        return self.ring.char()

    def order(self):
        return -1

    def setRepr(self, rep):
        self.repr = rep
        if rep == "reduced":
            self.ring.setRepr(rep)

    def __init__(self, ring: Ring, var, **kw):
            

        # For string conversion
        self.chain = 0
        if isinstance(ring,PolynomialRing):
            self.chain = ring.chain+1
        
        if var is None:
            self.var = VARS[self.chain]
        else:
            self.var = var

        self.ring = ring

        super().__init__(zero = self.Element([ring.zero]),one = self.Element([ring.one]))

        self.setRepr(self.ring.repr)
    
    def build(self,args):
        if isinstance(args[0], Ring.Element) and args[0].ring == self.ring:
            return self.Element(args)
        l = []
        for i in range(len(args)):
            l.append(self.ring.build(args[i]))
        return self.Element(l)
        


class PolynomialED(PolynomialRing, EuclideanDomain):
    """
        A polynomial ring which is an euclidean domain, or equivalently, whose ring of coefficients is a field
    """


    class Element(PolynomialRing.Element, EuclideanDomain.Element):

        def __floordiv__(self, other):
            super().__floordiv__(other)
            quot, rem = polynomial_division(self, other)
            return quot

        def __mod__(self, other):
            super().__mod__(other)
            quot, rem = polynomial_division(self, other)
            return rem

        def is_prime(self):
            if self.ring.ring.is_finite():
                return rabin_test(self)
            else:
                raise NotImplementedError()

        def is_unit(self):
            return self.deg() == 0 and self != self.ring.zero

        def factors(self):
            raise NotImplementedError()

        def croot(self):
            """ cth root of the element, where c is the characteristic of the ring """
            R = self.ring
            p = R.char()
            newc = []
            k = 0

            for i,c in enumerate(self.coefs):
                if i == k*p:
                    newc.append(c.croot())
                    k += 1

            return R.build(newc)

    def phi(self, element):
        return element.deg()

    def numphi(self, n: int):
        return self.ring.char() ** n



def polynomial_division(a, b):
    """ 
        Computes a // b, where a and b are polynomials.
        Returns (quotient, remainder)
        Long division algorithm
    """
    
    R = a.ring.ring

    if a.deg() < b.deg():
        return (a.ring.zero, a)

    # coef and exp of the largest monomials
    ca = a.deg()-1
    ea = a.coefs[-1]
    cb = b.deg()-1
    eb = b.coefs[-1]

    main_quot = a.__class__(([R.zero] * (ca-cb)) +  [ea / eb])
    reduced = a - b*main_quot

    quot, rem = polynomial_division(reduced, b)
    return (quot+main_quot, rem)


def factorize(n: int):
    # Crude factorization
    # TODO replace by sieve
    if n < 2:
        return []
    for i in range(2,n+1):
        if n % i == 0:
            return [i] + factorize(n//i)


def rabin_test(pol):
    """
        Rabin's test of irreducibility for polynomials with coefficients in a finite field
    """
    PR = pol.ring       # polinomial ring
    R = PR.ring         # underlying ring (field)
    assuming(R.is_finite())

    p = R.order()
    n = pol.deg()

    x = PR.build([R.zero,R.one])

    if (x**(p**n)-x) % pol != PR.zero:
        return False

    factors = set(factorize(n))
    for f in factors:
        ni = n//f
        if not gcd(pol, (x**(p**ni)-x)).is_unit():
            return False

    return True

def GetPolynomials(ring, var = None):
    if isinstance(ring, Field):
        return PolynomialED(ring, var)
    else:
        return PolynomialRing(ring, var)
