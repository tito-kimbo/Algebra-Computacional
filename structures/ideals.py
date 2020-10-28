from abc import ABC,abstractmethod
from structures.rings import Ring, EuclideanDomain
from structures.fields import Field
from algorithms.eea import eea
from utils import assuming

class Ideal(ABC):
    """Class representing an ideal over a ring."""
    def __init__(self,ring,generators):
        assuming(isinstance(ring, Ring),"ring must be a Ring")

        for g in generators:
            assuming(isinstance(g, ring.Element))

        self.ring = ring
        self.generators = generators
    
    @abstractmethod
    def has(self,element):
        pass
    
    @abstractmethod
    def is_principal(self):
        pass

    @abstractmethod
    def is_maximal(self):
        pass


class BaseQuotient(ABC):
    
    def __init__(self, ring, ideal, **kw):
        self.ring = ring
        self.ideal = ideal
        super().__init__(**kw)

    def build(self, *args, rep=None, **kwargs):
        if rep is not None:
            return self.Element(rep)
        else:
            return self.Element(self.ring.build(*args, **kwargs))


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
            
        def opp(self):
            return self.__class__(self.rep.opp())

        def reduce_rep(self):
            raise NotImplemented()
    
    def __init__(self,ring,ideal,**kw):

        assuming(isinstance(ring,Ring), "ring must be a Ring")
        assuming(isinstance(ideal,Ideal), "ideal must be an Ideal")
        
        super().__init__(ring=ring, ideal=ideal, zero=self.Element(ring.zero), one=self.Element(ring.one), **kw)


class FieldQuotient(BaseQuotient, Field):

    class Element(RingQuotient.Element, Field.Element):
        
        def inverse(self):
            raise NotImplemented()

        def __truediv__(self, other):
            super().__truediv__(self,other)
            return self * other.inverse

        def __mod__(self, other):
            super().__truediv__(self,other)
            return self - (self / other)


    def __init__(self,ring,ideal,inverse=None,**kw):

        assuming(ideal.is_maximal, "ideal must be maximal")
        
        if inverse:
            self.Element.inverse = inverse

        super().__init__(ring=ring, ideal=ideal, zero=self.Element(ring.zero), one=self.Element(ring.one), **kw)



# Existing quotients are stored to avoid duplication

_quotients = dict()

def _new_quotient(ring, ideal):
    if ideal.is_maximal():
        if isinstance(ring, EuclideanDomain):
            # We can calculate the inverse using the eea if the ring is an ED

            def inv(self):
                gcd, inv, _ = eea(self.rep,self.ring.ideal.generators[0])
                return self.__class__(inv)

            return FieldQuotient(ring=ring, ideal=ideal,inverse=inv)
        return FieldQuotient(ring, ideal)
    else:
        return RingQuotient(ring, ideal)

def Quotient(ring, ideal):
    if (ring,ideal) in _quotients:
        return _quotients[(_ring,ideal)]
    else:
        quot = _new_quotient(ring, ideal)
        _quotients[(ring,ideal)] = quot
        return quot
