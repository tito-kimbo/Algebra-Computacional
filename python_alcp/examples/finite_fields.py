from python_alcp.structures.ideals import GetQuotient
from python_alcp.structures.polynomials import GetPolynomials
from python_alcp.examples.rings import Z
from python_alcp.utils import assuming


def FiniteField(p, pol = None, var = None):
    if pol is None:
        return GetQuotient(Z, p*Z)
    else:
        Pols = GetPolynomials(Z/(p*Z), var)
        assuming(Pols(pol).is_prime(), f"{Pols(pol)} is not irreducible on {Pols}")
        return GetQuotient(Pols, Pols*Pols(pol))
