from abc import ABC,abstractmethod
from structures.rings import Ring, EuclideanDomain
from structures.fields import Field
from algorithms.divisibility import gcd, eea
from utils import assuming

class Ideal(ABC):
    """Class representing an ideal over a ring."""
    def __init__(self,ring,generators):
        assuming(isinstance(ring, Ring),"ring must be a Ring")

        for g in generators:
            assuming(isinstance(g, ring.Element))

        self.ring = ring
        self.generators = generators

    def __eq__(self,other):
        if isinstance(other, Ideal):
            # Better way?
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

    # Adds support for stuff like Z/NZ(5) instead of Quotient(Z, NZ(5))
    def __rtruediv__(self, other):
        return Quotient(other, self)

    def __str__(self):
        insides = ','.join(map(str, self.generators))
        return f"({insides})"


class EDIdeal(Ideal):
    """ Ideal of an EuclideanDomain """

    def __init__(self,ring, generators):
        assuming(isinstance(ring, EuclideanDomain),"ring must be an EuclideanDomain")

        for g in generators:
            assuming(isinstance(g, ring.Element))

        if len(generators) > 1:
            self.generator = gdc(generators)
        else:
            self.generator = generators[0]

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


class BaseQuotient(ABC):
    
    def __init__(self, ring, ideal, **kw):
        self.ring = ring
        self.ideal = ideal
        super().__init__(**kw)

    def __eq__(self,other):
        if isinstance(other, BaseQuotient):
            return self.ring == other.ring and self.ideal == other.ideal
        else:
            return False

    def build(self, *args, rep=None, **kwargs):
        if rep is not None:
            return self.Element(rep)
        else:
            return self.Element(self.ring.build(*args, **kwargs))

    def __str__(self):
        return f"{self.ring}/{self.ideal}"


class RingQuotient(BaseQuotient, Ring):

    class Element(Ring.Element): 
        """Class representing an element of the ring."""
        
        def __init__(self,rep):
            assuming(isinstance(rep,Ring.Element),
                    f"rep must be a ring element")

            self.rep = rep
        
        def __add__(self,other):
            super().__add__(other)
            return self.__class__(self.rep+other.rep)
            
        def __sub__(self,other):
            super().__sub__(other)
            return self.__class__(self.rep-other.rep)
        
        def __mul__(self,other):
            super().__mul__(other)
            return self.__class__(self.rep*other.rep)
        
        def __eq__(self,other):
            if not super().__eq__(other):
                return False
            return self.ring.ideal.has(self.rep-other.rep)
            
        def __str__(self):
            return "["+str(self.rep)+"]"
            
        def __neg__(self):
            return self.__class__(-self.rep)

        def reduce_rep(self):
            raise NotImplementedError()
    
    def __init__(self,ring,ideal,**kw):

        assuming(isinstance(ring,Ring), "ring must be a Ring")
        assuming(isinstance(ideal,Ideal), "ideal must be an Ideal")
        
        super().__init__(ring=ring, ideal=ideal, zero=self.Element(ring.zero), one=self.Element(ring.one), **kw)


class FieldQuotient(BaseQuotient, Field):

    class Element(RingQuotient.Element, Field.Element):
        
        def inverse(self):
            if self == self.ring.zero:
                raise ValueError(f"{self} does not have an inverse")
            return self._inverse()

        def _inverse(self):
            gcd, inv, _ = eea(self.rep,self.ring.ideal.generator)
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

        def is_prime(self):
            raise NotImplementedError()


    def __init__(self,ring,ideal,_inverse=None,**kw):

        assuming(ideal.is_maximal, "ideal must be maximal")
        
        if _inverse is not None:
            self.Element._inverse = _inverse

        super().__init__(ring=ring, ideal=ideal, zero=self.Element(ring.zero), one=self.Element(ring.one), **kw)



def Quotient(ring, ideal):
    if isinstance(ring, EuclideanDomain) and ideal.is_maximal():
        return FieldQuotient(ring, ideal)
    else:
        return RingQuotient(ring, ideal)
