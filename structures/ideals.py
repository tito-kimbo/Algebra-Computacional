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
class Quotient:

    def __init__(self,ring,ideal):
        assert isinstance(ring,Ring) and isinstance(ideal, Principal_Ideal)
        self.ring = ring
        self.ideal = ideal   
    