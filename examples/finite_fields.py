from structures.ideals import FieldQuotient
from structures.polynomials import GetPolynomials
from examples.rings import Z
from algorithms.discrete_log import discrete_log


class FiniteField(FieldQuotient):

    class Element(FieldQuotient.Element):

        def __str__(self):
            if self.ring.repr == "reduced":
                return super().__str__()
            elif self.ring.repr != None:
                if self == self.ring.zero:
                    return "0"
                if self == self.ring.one:
                    return "1"
                else:
                    dl = discrete_log(self.ring.generator(), self)
                    return f"{self.ring.repr}**{dl}"
            else:
                return super().__str__()
    
    def __init__(self, p, pol = None, var = None):
        if pol is None:
            super().__init__(Z, Z*p)
        else:
            Pols = GetPolynomials(Z/(p*Z), var)
            super().__init__(Pols, Pols*Pols.build(pol))
        self.setRepr("reduced")

    repr = None
    def setRepr(self, rep):
        self.repr = rep
        if rep == "reduced":
            self.ring.setRepr(rep)
