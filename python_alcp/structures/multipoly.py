from itertools import zip_longest

from python_alcp.structures.rings import (
        Ring,
        IntegralDomain,
        EuclideanDomain,
        RingElement,
        EuclideanDomainElement
)
from python_alcp.utils import (
        assuming,
        print_superscript,
        external,
        externals,
        prime_factors
)

class Monomial():
    def __init__(self, deg, vars):
        # TYPECHECKING
        self.deg = deg
        self.vars = vars
    
    def __repr__(self):
        return str(self)
    
    def __str__(self):
        return f''.join([''.join([str(y),'^',str(x)]) for x,y in zip(self.deg, self.vars) if x > 0])
        
    def __mul__(self,other):
        # We assume variables match!
        d = [self.deg[i]+other.deg[i] for i in range(len(self.deg))]
        return Monomial(d,self.vars)
        
    def __eq__(self,other):
        return self.deg == other.deg

class MultivariatePolynomial(RingElement):
    
    def __init__(self,coefs,monomials=[]):
        # TYPECHECKING IS IMPORTANT HERE
        if type(coefs) != list:
            coefs = [coefs]
        
        if not isinstance(coefs[0],self.coefRing):
            coefs = [self.coefRing(x) for x in coefs]
        self.coefs = coefs
        self.monomials = monomials
    
    def deg(self):
        return (max([x.deg[i] for x in self.monomials]) for i in range(len(self.monomials[0].deg)))
    
    def __add__(self,other):
        # WE NEED TO MATCH MONOMIALS
        return type(self)([x+y for x,y in zip_longest(self.val,other.val,fillvalue=type(self).coefRing.zero)])

    def __sub__(self,other):
        # WE NEED TO MATCH MONOMIALS
        if not hasattr(other, "val") or not hasattr(other.val, "__iter__"):
            other = type(self)(other)
        return type(self)([x-y for x,y in zip_longest(self.val,other.val,fillvalue=type(self).coefRing.zero)])
    
    # Convolutional product
    def inner_mul(self,other):
        # PRODUCT IS TOUGH 
        if not hasattr(other, "val") or not hasattr(other.val, "__iter__"):
            other = type(self)(other)
        R = type(self).coefRing
        newcs = [R.zero] * (self.deg() + other.deg() + 1)
        for i,v in enumerate(self.val):
            for j,w in enumerate(other.val):
                newcs[i+j] += v*w

        return type(self)(newcs)

    def __neg__(self):
        return MultivariatePolynomial([-x for x in self.coefs],self.monomials)

    def __floordiv__(self, other):
        # REALLY HARD
        if not hasattr(other, "val") or not hasattr(other.val, "__iter__"):
            other = type(self)(other)
        quot, rem = polynomial_division(self, other)
        return quot

    def __mod__(self, other):
        # REALLY HARD
        if not hasattr(other, "val") or not hasattr(other.val, "__iter__"):
            other = type(self)(other)
        quot, rem = polynomial_division(self, other)
        return rem

                
    def __eq__(self,other):
        # Monomials match and have same coefficients
        return hasattr(other, "val") and self.val == other.val
    
    def __str__(self):
        # TBD
        return self.coefs+self.monomials


    

    def inverse(self):
        assuming(self.is_unit(), "Can't invert non-unit")
        return type(self)(self.val[0].inverse())

    def der(self):
        # derivative
        new_val = []
        for i,c in enumerate(self.val[1:]):
            if c != type(self).zero:
                new_val.append((i+1)*c)
        return type(self)(new_val)


    def is_unit(self):
        return self.deg() == 0 and self.coefs[0].is_unit()

    def normal(self):
        raise NotImplementedError()
    
    def content(self):
        raise NotImplementedError()

    def primitive_part(self):
        raise NotImplementedError()

    def is_prime(self):
        raise NotImplementedError()

    def __hash__(self):
        return hash((type(self).__name__, self.coefs, self.monomials))

# Polynomials over integral domains
class MultivariatePolynomialRing(Ring):
    
    def __eq__(cls, other):
        return hasattr(other, "coefRing") and cls.coefRing == other.coefRing and ...
        hasattr(other, "vars") and cls.vars == other.vars

    def char(cls):
        return cls.coefRing.char()

    def order(cls):
        return 0

    def setRepr(cls, rep):
        cls.repr = rep
        if rep == "reduced":
            cls.coefRing.setRepr(rep)

    def generators(cls):
        # Could return all monomials of degree 1 such that the total sum of the multidegree is 1 + ring generators
        raise NotImplementedError()

    def units(cls):
        return {cls(u) for u in cls.coefRing.units()}

    def is_polynomial(cls):
        return True
    
    def monomial(cls,deg):
        return Monomial(deg,cls.vars)
        
def GetMultiPoly(ring, vars):
    attrs = {"coefRing": ring, "vars": vars}
    return MultivariatePolynomialRing(f"{ring}[{vars}]", (MultivariatePolynomial,), attrs)