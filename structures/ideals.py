from abc import ABC,abstractmethod
from rings import * 

# WIP
class Principal_Ideal(ABC):
	
    def __init__(self,ring,generator):
        self.ring = ring
        assert isinstance(self.generator, Ring.Element)
        self.generator = generator
    
    @abstractmethod
    def has(self,element):
        pass

# WIP
class Quotient(Ring):
    
    class QuotientElement(Ring.Element): 
        """Class representing an element of the ring."""
        
        def __init__(self,rep):
            assert isinstance(rep,Ring.Element)
            self.rep = rep
        
        def __add__(self,other):
            return QuotientElement(self.rep+other.rep)
            
        def __sub__(self,other):
            return QuotientElement(self.rep-other.rep)
        
        def __mul__(self,other):
            return QuotientElement(self.rep*other.rep)
        
        @abstractmethod
        def __eq__(self,other):
            pass # Usually it should just be checking for divisibility in the case of a principal ideal
            
        def __str__(self):
            return "["+str(self.rep)+"]"
            
        def opp(self):
            return QuotientElement(self.rep.opp())
    
    def __init__(self,ring,generator):
        assert isinstance(ring,Ring) and isinstance(generator, Ring.Element)
        self.ring = ring
        self.ideal = Principal_Ideal(ring,generator)
    