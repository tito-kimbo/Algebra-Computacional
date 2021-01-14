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

"""
Returns True if m1 divides m2. False otherwise. Assumes the
monomials are expressed over the same variables.

:m1 first monomial
:m2 second monomial
"""
def divides_monomial(m1,m2):
    return all(m1.deg[i] < m2.deg[i]  for i in range(len(m1.deg)))

class Monomial():
    def __init__(self, deg, vars):
        # TYPECHECKING    
        self.deg = tuple(deg)
        self.vars = vars
    
    def __repr__(self):
        return str(self)
    
    def __str__(self):
        return f''.join([''.join([str(y),'^',str(x)]) for x,y in zip(self.deg, self.vars) if x > 0])
        
    def __mul__(self,other):
        # We assume variables match!
        d = [self.deg[i]+other.deg[i] for i in range(len(self.deg))]
        return Monomial(d,self.vars)
    
    def __truediv__(self,other):
        pass
        
    def __eq__(self,other):
        return self.deg == other.deg
    
    def __hash__(self):
        return hash((type(self).__name__, self.deg))
    
    def is_zero(self):
        return all(x == 0 for x in self.deg)

class MultivariatePolynomial(RingElement):
    
    def __init__(self,coefs,monomials=[],is_dict=False):
        # TYPECHECKING IS IMPORTANT HERE
        if is_dict:
            self.coefs = coefs
        else:
            if type(coefs) != list:
                coefs = [coefs]
            if monomials == []:
                v = self.vars
                monomials = [Monomial((0 for _ in range(len(v))),v)]# WE NEED TO CREATE THE "ZERO" MONOMIAL
            
            if not isinstance(coefs[0],self.coefRing):
                coefs = [self.coefRing(x) for x in coefs]
            self.coefs = dict(zip(monomials,coefs))
    
    def deg(self):
        monomials = self.coefs.keys()
        return (max([x.deg[i] for x in monomials]) for i in range(len(monomials[0].deg)))
    
    def __add__(self,other):
        c = dict(self.coefs)
        for k in other.coefs:
            if k in c:
                c[k] += other.coefs[k]
            else:
                c[k] = other.coefs[k]
        return type(self)(c,is_dict=True)

    def __sub__(self,other):
        c = dict(self.coefs)
        for k in other.coefs:
            if k in c:
                c[k] -= other.coefs[k]
            else:
                c[k] = -other.coefs[k]
        return type(self)(c,is_dict=True)
    
    # Convolutional product
    def inner_mul(self,other):
        c = dict()
        for k1,v1 in self.coefs.items():
            for k2,v2 in other.coefs.items():
                k = k1*k2
                v = v1*v2
                if k in c:
                    c[k] += v
                else:
                    c[k] = v
        return type(self)(c,is_dict=True)

    def __neg__(self):
        return type(self)([-x for x in self.coefs],self.monomials)

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
        return hasattr(other, "val") and self.coefs == other.coefs
    
    def __str__(self):
        s = ''
        for k,v in self.coefs.items():
            if k.is_zero():
                s = s + ''.join([str(v), ' + '])
            else:
                s = s + ''.join([str(v),str(k), ' + '])
        s = s[:-3:]
        return s
    
    def __repr__(self):
        return self.__str__()
        
    def inverse(self):
        # This does not work
        assuming(self.is_unit(), "Can't invert non-unit")
        return type(self)(self.coefs[0].inverse())
    
    def is_unit(self):
        # This does not work
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