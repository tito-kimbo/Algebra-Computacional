from structures.ideals import FieldQuotient
from structures.polynomials import GetPolynomials
from examples.rings import Z


class FiniteField(FieldQuotient):

    
    def __init__(self, p, generator):
        Pols = GetPolynomials(Z/(p*Z))

        super().__init__(Pols, Pols*Pols.build(generator))

    
