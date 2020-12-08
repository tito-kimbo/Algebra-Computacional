from structures.ideals import FieldQuotient
from structures.polynomials import GetPolynomials
from examples.rings import Z


class FiniteField(FieldQuotient):

    
    def __init__(self, p, pol):
        Pols = GetPolynomials(Z/(p*Z))

        super().__init__(Pols, Pols*Pols.build(pol))

    def generator(self):
        """
            Returns the class x + <pol>, which is a generator of the field if pol is primitive
        """

        return self.build([0,1])
