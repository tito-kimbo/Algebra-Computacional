from abc import ABC,abstractmethod

class Ring(ABC):
    """Class representing an algebraic ring."""
    
    class Element(ABC): 
        """Class representing an element of the ring."""
        @abstractmethod
        def __add__(self,other):
            pass
        @abstractmethod
        def __sub__(self,other):
            pass
        @abstractmethod
        def __mul__(self,other):
            pass
        @abstractmethod
        def __eq__(self,other):
            pass
        @abstractmethod
        def __str__(self):
            pass
        @abstractmethod
        def opp(self):
            pass

    def __init__(self,zero,one,elementClass=None):
        self.zero = zero
        self.one = one
        if elementClass is None:
            self.elementClass = Ring.Element
        else:
            self.elementClass = elementClass

class IntegralDomain(Ring):
    """Class representing an integral domain."""

    class IDElement(Ring.Element): 
        """Class representing an element of the ID."""
        @abstractmethod
        def __truediv__(self,other):
            pass
        @abstractmethod
        def __mod__(self,other):
            pass
    
    def __init__(self,zero,one,elementClass=None):
        if elementClass is None:
            elementClass = IntegralDomain.IDElement
        super().__init__(zero,one,elementClass)
    
class EuclideanDomain(IntegralDomain):
    """Structure representing an euclidean domain."""
    
    def __init__(self,zero,one,elementClass=None):
        if elementClass is None:
            elementClass = IntegralDomain.IDElement
        super().__init__(zero,one,elementClass)
    
    @abstractmethod
    def phi(self,element):
        pass