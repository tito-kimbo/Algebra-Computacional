from abc import ABC,abstractmethod
from utils import assuming

class Ring(ABC):
    """Class representing an algebraic ring (commutative ring with unity)"""
    
    class Element(ABC): 
        """Class representing an element of the ring."""

        @abstractmethod
        def __add__(self,other):
            assuming(self.ring == other.ring,
                    "You can only add elements of the same ring")
            pass

        @abstractmethod
        def __sub__(self,other):
            assuming(self.ring == other.ring,
                    "You can only substract elements of the same ring")
            pass

        @abstractmethod
        def __mul__(self,other):
            assuming(self.ring == other.ring,
                    "You can only multiply elements of the same ring")
            pass

        @abstractmethod
        def __eq__(self,other):
            return self.ring == other.ring
            
        @abstractmethod
        def __str__(self):
            pass

        @abstractmethod
        def opp(self):
            pass


    def __init__(self,zero,one, **kw):
        self.zero = zero
        self.one = one
        self.Element.ring = self
        super().__init__(**kw)
    
    def build(self,*args,**kwargs):
        return self.Element(*args,**kwargs)

    @abstractmethod
    def __eq__(self, other):
        pass

    @abstractmethod
    def __str__(self):
        pass


class IntegralDomain(Ring):
    """Class representing an integral domain."""

    class Element(Ring.Element): 
        """Class representing an element of the ID."""
        @abstractmethod
        def __truediv__(self,other):
            assuming(self.ring == other.ring,
                    "You can only divide elements of the same ring")
            pass
        @abstractmethod
        def __mod__(self,other):
            assuming(self.ring == other.ring,
                    "You can only calculate the remainder between elements of the same ring")
            pass
    
    
class EuclideanDomain(IntegralDomain):
    """Structure representing an euclidean domain."""
    
    @abstractmethod
    def phi(self,element):
        pass
