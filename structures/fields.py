from abc import ABC,abstractmethod
from structures.rings import EuclideanDomain, IntegralDomain

class Field(EuclideanDomain):
    """Class representing a field."""
    
    class Element(IntegralDomain.Element):
        """Class representing an element of a field."""
        @abstractmethod
        def inverse(self):
            pass

        @abstractmethod
        def __truediv__(self,other):
            pass
        
    def phi(self,element):
        return 1;
