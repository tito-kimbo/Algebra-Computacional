from python_alcp.structures.rings import EuclideanDomain, IntegralDomain, UniqueFactorizationDomain
from python_alcp.structures.polynomials import PolynomialRing, polynomial_division
from python_alcp.algorithms.divisibility import gcd

def gcd_dfu(f, g, *args):
    """
        f and g must be polynomials with coefficients in a UFD
    """
    
    # We check both elements belong to the same ring
    if f.ring != g.ring:
        raise ValueError("f and g must be elements of the same ring")

    # And the coefficient ring is an UFD
    if not f.ring.coefRing.is_ufd():
        raise ValueError("The elements belong to a ring that is not an UFD")
    
    # We put both arguments in their normal form
    f = f.normal()
    g = g.normal()

    # First, if deg(f) < deg(g), we have to swap both elements
    if f.deg() < g.deg():
        f, g = g, f

    R =  f.ring
    r = [f,g]

    while r[-1] != R.zero:
        # Then, we take the remainder of the division
        quot, re = polynomial_division(r[-2], r[-1], pseudo = True)
        
        # Lastly, we take the primitive part of the polynomial
        re = re.primitive_part()

        r.append(re)

    # If we are calculating the gcd of two elements, we process the next element
    if len(args) > 0:
        return gcd_dfu(r[-2], *args)
    else:
    #If not, we just return the normal form of the result
        return r[-2].normal()
