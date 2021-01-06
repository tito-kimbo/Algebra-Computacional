from functools import reduce
from abc import ABC, abstractmethod
import random

from python_alcp.utils import assuming
from python_alcp.structures.rings import Ring, RingElement, FieldElement, Field
from python_alcp.utils import external, externals

class Ideal(ABC):
    """Class representing an ideal over a ring."""

    def __init__(self, ring, generators):
        self.ring = ring
        self.generators = [ring(g) for g in generators]

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

    def __init__(self, ring, generators):

        assuming(len(generators) > 0, "At least one generator is needed")
        assuming(ring.is_euclidean(), "Ring must be Euclidean")

        gen = [ring(g) for g in generators]

        # We can calculate the generator as the gcd of all the given elements

        if len(generators) > 1:
            self.generator = externals.gcd(*gen)
        else:
            self.generator = gen[0]

        super().__init__(ring, [self.generator])


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


@external
def GetIdeal(ring, generators):
    if ring.is_euclidean():
        return EDIdeal(ring, generators)
    else:
        return Ideal(ring, generators)



class Quotient(Ring):

    """ Abstract superclass for quotients """

    def __eq__(cls,other):
        return (hasattr(other, "baseRing") and hasattr(other, "ideal")
                and cls.baseRing == other.baseRing and cls.ideal == other.ideal)

    def __str__(cls):
        return f"({cls.baseRing}/{cls.ideal})"

    def char(cls):
        if hasattr(cls, "_char"):
            return cls._char

        c = 1
        acc = cls.one
        while c < 1000 * 1000:
            if acc == cls.zero:
                cls._char = c
                return c
            acc = acc + cls.one
            c += 1

        return 0

    def order(cls):
        R = cls.baseRing
        n = R.numphi(R.phi(cls.ideal.generator))
        if n > 0:
            return n
        else:
            return 0

    def is_finite(cls):
        return cls.order() != 0

    def generators(cls):
        return {cls(g) for g in cls.baseRing.generators() if cls(g) != cls.zero}

    def units(cls):
        return {e for e in cls.elements() if e.is_unit()}

    def generator(cls):
        """
            Returns a generator of the field, if one can be found.
            Returns [x] in the case of finite fields
        """

        if cls.baseRing.is_polynomial():
            return cls.build([0,1])
        elif cls.char() > 2:
            return cls.build(2)
        else:
            return cls.one

class FieldQuotient(Quotient, Field):
    pass

class RingQuotientElement(RingElement):

    """ A quotient which has a ring structure """

    def __init__(self, val):
        return super().__init__(self.reduce_rep(self.baseRing(val)))
    
    def __add__(self,other):
        return type(self)(self.val+other.val)
        
    def __sub__(self,other):
        return type(self)(self.val-other.val)
    
    def inner_mul(self,other):
        return type(self)(self.val*other.val)

    def __floordiv__(self, other):
        return type(self)(self.val // other.val)

    def __truediv__(self,other):
        return self * other.inverse()

    def __mod__(self, other):
        return type(self)(self.val % other.val)
    
    def __eq__(self,other):
        return type(other) == type(self) and self.ideal.has(self.val-other.val)

    def inverse(self):
        return type(self)(externals.modinv(self.val, type(self).ideal.generator))
        
    def __str__(self):
        if self.ring.repr == "reduced":
            return str(self.val)
        else:
            return f"[{self.val}]"
        
    def __neg__(self):
        return type(self)(-self.val)

    @classmethod
    def reduce_rep(cls, rep):
        if hasattr(rep, "__mod__"):
            return rep % cls.ideal.generator
        else:
            return rep

    def is_unit(self):
        return externals.gcd(self.ring.ideal.generator, self.val).is_unit()

    def __lt__(self, other):
        return self.reduce_rep(self.val) < other.reduce_rep(other.val)

    def __hash__(self):
        return hash((type(self).__name__, self.val))


class FieldQuotientElement(RingQuotientElement, FieldElement):
    """ A quotient which has a field structure """

    def __truediv__(self, other):
        return self * other.inverse()

    def __floordiv__(self, other):
        return self.__truediv__(other)

    def __mod__(self, other):
        return self.ring.zero

    def croot(self):
        """
            Computes the cth root of the element, where c is the characteristic
            of the ring
        """
        if self == self.ring.zero:
            return self

        # Assumes generator is primitive
        # Multiplicative ring order
        if self.ring.char() == self.ring.order():
            return self

        Z = externals.Z

        o = self.ring.order() - 1
        p = self.ring.char()
        g = self.ring.build([0,1])
        log = externals.discrete_log(g, self)
        assuming(log is not None, "The chosen field generator is not primitive")
        # Find an integer x such that x*p = log (mod o)
        gcd, __, inv = externals.eea(Z.build(o),Z.build(p))
        if gcd != Z.one:
            assuming(gcd.is_unit(), f"Could not find c-th root of {self} in {type(self)}. Is the field well constructed?")
            inv = inv // gcd
        x = inv*Z.build(log) % Z.build(o)
        return g**(x.val)

    def is_prime(self):
        raise NotImplementedError()

@external
def GetQuotient(ring, ideal):
    if ring.is_euclidean() and ideal.is_maximal():
        return GetFieldQuotient(ring, ideal)
    else:
        return GetRingQuotient(ring, ideal)

@external
def GetRingQuotient(ring, ideal):
    bases = (RingQuotientElement,)
    return Quotient(f"{ring}/{ideal}", bases, {'baseRing': ring, 'ideal': ideal})


def GetFieldQuotient(ring, ideal):
    bases = (FieldQuotientElement,)
    return FieldQuotient(f"{ring}/{ideal}", bases, {'baseRing': ring, 'ideal': ideal})

