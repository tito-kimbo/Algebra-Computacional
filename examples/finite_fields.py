from structures.ideals import FieldQuotient
from structures.polynomials import GetPolynomials
from examples.rings import Z
from algorithms.discrete_log import discrete_log


class FiniteField(FieldQuotient):

    repr = None

    class Element(FieldQuotient.Element):

        def __str__(self):
            if self.ring.repr != None:
                if self == self.ring.zero:
                    return "0"
                if self == self.ring.one:
                    return "1"
                else:
                    dl = discrete_log(self.ring.generator(), self)
                    return f"{self.ring.repr}**{dl}"
            else:
                return super().__str__()
    
    def __init__(self, p, pol):
        Pols = GetPolynomials(Z/(p*Z))

        super().__init__(Pols, Pols*Pols.build(pol))

    def setRepr(self, rep):
        self.repr = rep

    def generator(self):
        """
            Returns the class x + <pol>, which is a generator of the field if pol is primitive
        """

        return self.build([0,1])
