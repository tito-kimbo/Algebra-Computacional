from structures.ideals import EDIdeal, FieldQuotient
from structures.polynomials import PolynomialRing


# WIP

class FiniteField(FieldQuotient):

    class Element(FieldQuotient.Element):

        def derivative(self):
    
    
    def __init__(self, p, generator):
        super().__init__(PolynomialRing(Z/(p*Z)), EDIdeal(generator.ring, [generator]))

    
