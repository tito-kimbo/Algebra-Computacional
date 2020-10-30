from abc import ABC,abstractmethod
from utils import assuming

class Ring(ABC):
    """Class representing an algebraic ring (commutative ring with unity)"""
    
    class Element(ABC): 
        """Class representing an element of the ring."""

        ring = None     # See Ring.__new__

        @abstractmethod
        def __add__(self,other):
            assuming(isinstance(other, self.ring.Element) and self.ring == other.ring,
                    "You can only add elements of the same ring")
            pass

        @abstractmethod
        def __sub__(self,other):
            assuming(isinstance(other, self.ring.Element) and self.ring == other.ring,
                    "You can only substract elements of the same ring")
            pass

        @abstractmethod
        def __mul__(self,other):
            assuming(isinstance(other, self.ring.Element) and self.ring == other.ring,
                    "You can only multiply elements of the same ring")
            pass

        @abstractmethod
        def __eq__(self,other):
            return isinstance(other, self.__class__) and self.ring == other.ring
            
        @abstractmethod
        def __str__(self):
            pass
        
        def __repr__(self):
            return self.__str__()

        @abstractmethod
        def __neg__(self):
            pass

        def __repr__(self):
            return self.__str__()


    def __new__(cls, *args, **kwargs): 
        # Decorates the Element class before construction of the ring
        # so that self.Element.ring == self

        obj = super(Ring, cls).__new__(cls)

        obj.Element = type('.'.join([cls.__module__,cls.Element.__qualname__]), (cls.Element,), {"ring": obj})
        return obj

    def __init__(self,zero,one, **kw):
        self.zero = zero
        self.one = one

        super().__init__(**kw)
    
    def build(self,*args,**kwargs):
        return self.Element(*args,**kwargs)

    @abstractmethod
    def __eq__(self, other):
        pass
    
    @abstractmethod
    def __str__(self):
        pass
    
    def __repr__(self):
        return self.__str__()


class IntegralDomain(Ring):
    """Class representing an integral domain."""

    class Element(Ring.Element): 
        """Class representing an element of the ID."""
        @abstractmethod
        def __floordiv__(self,other):
            assuming(self.ring == other.ring,
                    "You can only divide elements of the same ring")
            pass
        @abstractmethod
        def __mod__(self,other):
            assuming(self.ring == other.ring,
                    "You can only calculate the remainder between elements of the same ring")
            pass

        @abstractmethod
        def is_prime(self):
            pass
    
    
class EuclideanDomain(IntegralDomain):
    """Structure representing an euclidean domain."""
    
    @abstractmethod
    def phi(self,element):
        pass
