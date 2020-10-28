from abc import ABC,abstractmethod
from rings import EuclideanDomain, IntegralDomain

class Field(EuclideanDomain):
    """Class representing a field."""
    
    class FieldElement(IntegralDomain.IDElement):
        """Class representing an element of a field."""
        @abstractmethod
        def inverse(self):
            pass
        
    def phi(self,element):
        return 1;
    
    def __init__(self,zero,one,elementClass=None):
        if elementClass is None:
            elementClass = Field.FieldElement
        super().__init__(zero,one,elementClass)


