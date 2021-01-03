from abc import ABC,abstractmethod
from structures.rings import Ring
from structures.fields import Field
from algorithms.divisibility import gcd, eea
from utils import assuming
from functools import reduce

class Ideal(ABC):
    """Class representing an ideal over a ring."""

    def __init__(self,generators: Ring.Element):
        self.ring = generators[0].ring
        self.generators = generators

    def __eq__(self,other):
        if isinstance(other, Ideal):
            return self.ring == self.ring and self.generators == other.generators
        else:
            return False
    
    @abstractmethod
    def has(self,element):
        pass
    
    @abstractmethod
    def is_principal(self):
        pass

    @abstractmethod
    def is_maximal(self):
        pass

    def __str__(self):
        insides = ','.join(map(str, self.generators))
        return f"({insides})"

    def __repr__(self):
        return self.__str__()


class EDIdeal(Ideal):
    """ Ideal of an EuclideanDomain """

    def __init__(self, generators):

        assuming(len(generators) > 0, "At least one generator is needed")
        assuming(generators[0].ring.is_euclidean(), "Ring must be Euclidean")

        # We can calculate the generator as the gcd of all the given elements

        if len(generators) > 1:
            self.generator = gcd(generators)
        else:
            self.generator = generators[0]

        super().__init__([self.generator])


    def has(self,element):
        return element % self.generator == self.ring.zero

    def is_principal(self):
        # Euclidean domains are principal ideal domains
        return True

    def is_maximal(self):
        # In a PID, maximal <=> nonzero and prime
        # Also, nonzero and prime <=> ideal = (p) and p is prime
        return self.generator.is_prime()

    def __str__(self):
        return f"{self.generator}{self.ring}"


def GetIdeal(generators):
    g = generators[0]
    if g.ring.is_euclidean():
        return EDIdeal(generators)
    else:
        return Ideal(generators)



class BaseQuotient(ABC):

    """ Abstract superclass for quotients """

    def __init__(self, ring: Ring, ideal: Ideal, **kw):
        self.ring = ring
        self.ideal = ideal
        super().__init__(**kw)

    def __eq__(self,other):
        if isinstance(other, BaseQuotient):
            return self.ring == other.ring and self.ideal == other.ideal
        else:
            return False

    def build(self, *args, rep: Ring.Element = None, **kwargs):
        if rep is not None:
            return super().build(rep)
        else:
            return super().build(self.ring.build(*args, **kwargs))

    def __str__(self):
        return f"{self.ring}/{self.ideal}"


    def char(self):
        if "_char" in self.__dict__:
            return self._char

        c = 1
        acc = self.one
        while c < 1000 * 1000:
            if acc == self.zero:
                self._char = c
                return c
            acc = acc + self.one
            c += 1

    def order(self):
        R = self.ring
        n = R.numphi(R.phi(self.ideal.generator))
        if n > 0:
            return n
        else:
            return 0


class RingQuotient(BaseQuotient, Ring):

    """ A quotient which has a ring structure """

    class Element(Ring.Element): 
        """Class representing an element of the ring."""
        
        def __init__(self, rep: Ring.Element):
            assuming(rep.ring == self.ring.ring,
                    f"rep must be a ring element")
            if rep.ring.is_euclidean():
                rep = self.reduce_rep(rep)
            self.val = rep
        
        def __add__(self,other):
            super().__add__(other)
            return self.__class__(self.val+other.val)
            
        def __sub__(self,other):
            super().__sub__(other)
            return self.__class__(self.val-other.val)
        
        def inner_mul(self,other):
            return self.__class__(self.val*other.val)
        
        def __eq__(self,other):
            return isinstance(other, Ring.Element) and other.ring == self.ring and self.ring.ideal.has(self.val-other.val)
            
        def __str__(self):
            if self.ring.repr == "reduced":
                return str(self.val)
            else:
                return f"[{self.val}]"
            
        def __neg__(self):
            return self.__class__(-self.val)

        def reduce_rep(self, rep):
            assuming(rep.ring.is_euclidean(), "Can't compute reduced rep in non-euclidean domain")
            return rep % self.ring.ideal.generator

        def is_unit(self):
            raise NotImplementedError()

        def __lt__(self, other):
            return self.reduce_rep(self.val) < other.reduce_rep(self.val)
    
    repr = None

    def setRepr(self, rep):
        self.repr = rep

    def __init__(self, ring: Ring, ideal: Ideal, **kw):
        self.ring = ring
        self.ideal = ideal
        super().__init__(ring=ring, ideal=ideal, zero=self.Element(ring.zero), one=self.Element(ring.one), **kw)



class FieldQuotient(BaseQuotient, Field):
    """ A quotient which has a field structure """

    class Element(RingQuotient.Element, Field.Element):

        def __init__(self, rep: Ring.Element):
            if rep.ring.is_euclidean():
                rep = self.reduce_rep(rep)
            super().__init__(rep=rep)
        
        def inverse(self):
            assuming(self.ring.is_euclidean(), "Can't calculate modular inverse in non-euclidean domain")
            if self == self.ring.zero:
                raise ValueError(f"{self} does not have an inverse")
            gcd, inv, _ = eea(self.val,self.ring.ideal.generator)
            return self.__class__(inv)


        def __truediv__(self, other):
            super().__truediv__(other)
            return self * other.inverse()

        def __floordiv__(self, other):
            super().__floordiv__(other)
            return self.__truediv__(other)

        def __mod__(self, other):
            super().__mod__(other)
            return self.ring.zero

        def reduce_rep(self, rep):
            assuming(rep.ring.is_euclidean(), "Can't compute reduced rep in non-euclidean domain")
            return rep % self.ring.ideal.generator

        def croot(self):
            """
                Computes the cth root of the element, where c is the characteristic
                of the ring
            """

            if self == self.ring.zero:
                return self

            from algorithms.divisibility import eea
            from algorithms.discrete_log import discrete_log
            from examples.rings import Z
            # Assumes generator is primitive
            # Multiplicative ring order
            if self.ring.char() == self.ring.order():
                return self

            o = self.ring.order() - 1
            p = self.ring.char()
            g = self.ring.build([0,1])
            log = discrete_log(g, self)
            assuming(log is not None, "The chosen field generator is not primitive")
            # Find an integer x such that x*p = log (mod o)
            _, __, inv = eea(Z.build(o),Z.build(p))
            x = inv*Z.build(log) % Z.build(o)
            return g**(x.val)


        def is_prime(self):
            raise NotImplementedError()


    def is_finite(self):
        # TODO: This is not true in general
        return True


    def __init__(self,ring: Ring,ideal: Ideal,**kw):

        assuming(ideal.is_maximal(), "ideal must be maximal")

        self.ring = ring
        self.ideal = ideal
        super(Field, self).__init__(zero=self.Element(ring.zero), one=self.Element(ring.one), **kw)



def GetQuotient(ring, ideal):
    if ring.is_euclidean() and ideal.is_maximal():
        return FieldQuotient(ring, ideal)
    else:
        return RingQuotient(ring, ideal)
