from abc import ABC,abstractmethod
from structures.rings import Ring
from utils import assuming

class Ideal(ABC):
    """Class representing an ideal over a ring."""
    def __init__(self,ring,generators):
        assuming(isinstance(ring, Ring),"ring must be a Ring")

        #for g in generators:
        #    assert isinstance(g, ring.elementClass)
        self.ring = ring
        self.generators = generators
    
    @abstractmethod
    def has(self,element):
        pass
    
    def is_principal(self):
        return (isinstance(self.generators,list) and len(self.generators)) == 1

class Quotient(Ring):
    
    def __init__(self,ring,ideal):
        assuming(isinstance(ring,Ring), "ring must be a Ring")
        assuming(isinstance(ideal,Ideal), "ideal must be an Ideal")

        self.ring = ring
        self.ideal = ideal
        
        class QuotientElement(Ring.Element): 
            """Class representing an element of the ring."""
            
            def __init__(self,rep):
                assuming(isinstance(rep,Ring.Element), "rep must be a ring element")
                self.rep = rep
            
            def __add__(self,other):
                return QuotientElement(self.rep+other.rep)
                
            def __sub__(self,other):
                return QuotientElement(self.rep-other.rep)
            
            def __mul__(self,other):
                return QuotientElement(self.rep*other.rep)
            
            def __eq__(self,other):
                return ideal.has(self.rep-other.rep)
                
            def __str__(self):
                return "["+str(self.rep)+"]"
                
            def opp(self):
                return QuotientElement(self.rep.opp())
            
        self.zero = QuotientElement(ring.zero)
        self.one = QuotientElement(ring.one)
        self.elementClass = QuotientElement
    
    def build(self,*args,**kwargs):
        rep = self.ring.build(*args,**kwargs)
        return self.elementClass(rep)
