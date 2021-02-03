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

:m1 monomial
:m2 monomial
"""
def divides_monomial(m1,m2):
    return all(m1.deg[i] <= m2.deg[i]  for i in range(len(m1.deg)))

"""
Returns a tuple with the differences of the degrees of two monomials.

:m1 monomial
:m2 monomial
"""
def deg_diff(m1,m2):
    return tuple([m1.deg[i]-m2.deg[i] for i in range(len(m1.deg))])


"""
Returns the total degree of a mononial.

:m monomial
"""
def total_deg(m):
    return sum(m.deg)

########### MONOMIAL ORDERS ##############

"""
Less than operation for the lexicographic monomial order.
"""
def lt_lex(m1,m2):
    dd = deg_diff(m1,m2)
    i=0
    while i < len(dd) and dd[i] == 0:
        i += 1
    return i < len(dd) and dd[i]<0

"""
Less than operation for the reverse lexicographic monomial order.
"""
def lt_reverse_lex(m1,m2):
    dd = deg_diff(m1,m2)
    i=len(dd)-1
    while i>=0 and dd[i]==0:
        i -= 1
    return i>=0 and dd[i]>0

"""
Less than operation for the graded lexicographic monomial order.
"""
def lt_graded_lex(m1,m2):
    d1,d2 = total_deg(m1),total_deg(m2)
    return d1 < d2 or (d1 == d2 and lt_lex(m1,m2)) 

"""
Less than operation for the graded reverse lexicographic monomial order.
"""
def lt_graded_reverse_lex(m1,m2):
    d1,d2 = total_deg(m1),total_deg(m2)
    return d1 < d2 or (d1 == d2 and lt_reverse_lex(m1,m2))
    
##########################################

"""
Class representing a Monomial over a set of variables.
"""
class Monomial():
    def __init__(self, deg, vars):
        # TYPECHECKING    
        self.deg = tuple(deg)
        self.vars = vars
    
    def __repr__(self):
        return str(self)
    
    def __str__(self):
        return ''.join([''.join([str(y),print_superscript(x)]) for x,y in zip(self.deg, self.vars) if x > 0])
        
    def __mul__(self,other):
        # We assume variables match!
        d = [self.deg[i]+other.deg[i] for i in range(len(self.deg))]
        return Monomial(d,self.vars)
    
    def __truediv__(self,other):
        if divides_monomial(other,self):
            return Monomial(deg_diff(self,other),self.vars)
        else:
            raise ValueError('Invalid division.')
        
    def __eq__(self,other):
        return self.deg == other.deg
       
    def __lt__(self,other):
        return lt_graded_lex(self,other)
    
    def __hash__(self):
        return hash((type(self).__name__, self.deg))
    
    def is_zero(self):
        return all(x == 0 for x in self.deg)

"""
Utility functon that divides two tuples element-by-element.
"""
def tuple_div(t1,t2):
    return tuple([t1[i]/t2[i] for i in range(len(t1))])

"""
Utility functio to turn a term of a polynomial into a 
"""
def term_to_poly(t,p_type):
    return p_type({t[0] : t[1]},is_dict=True)

"""
Chooses a polynomial from F such that its leading term divides the leading term of p.

:p  polynomial
:F  list of polynomials in the same polynomial ring as p
"""
def _choose_lt_divisor(p,F):
    p_lt = p.lt()
    f_lt = [f.lt() for f in F]
    idx = 0
    while idx<len(F) and not divides_monomial(f_lt[idx][0], p_lt[0]):
        idx += 1
    return idx

"""
Returns the quotient and modulo of p/F in reduced normal form.

:p dividend polynomial
:F polynomial or list of polynomials to use as divisors
"""
def div_poly_RNF(p,F):
    if type(F) != list:
        F = [F]
    h=type(p)(p.coefs,is_dict=True)
    A=[type(p)(0) for _ in F]
    r=type(p)(0)
    while h != type(p).zero:
        l = _choose_lt_divisor(h,F)
        while l < len(F):
            t = tuple_div(h.lt(),F[l].lt())
            aux_poly = term_to_poly(t,type(p))
            A[l] = A[l] + aux_poly
            h = h - aux_poly*F[l]
            l = _choose_lt_divisor(h,F)
        aux_poly = term_to_poly(h.lt(),type(p))
        r=r+aux_poly
        h =h-aux_poly
    return (A,r)

"""
Class representing a multivariate polynomial.
"""
class MultivariatePolynomial(RingElement):
    
    def __init__(self,coefs,monomials=[],is_dict=False):
        # TYPECHECKING IS IMPORTANT HERE
        if not is_dict:
            if type(coefs) != list:
                coefs = [coefs]
            if monomials == []:
                v = self.vars
                monomials = [Monomial((0 for _ in range(len(v))),v)]# WE NEED TO CREATE THE "ZERO" MONOMIAL
            
            if not isinstance(coefs[0],self.coefRing):
                coefs = [self.coefRing(x) for x in coefs]
            coefs = dict(zip(monomials,coefs))    
        
        # Reduce the representation by removing 0s
        self.coefs = {k : coefs[k] for k in sorted(coefs,reverse=True) if coefs[k] != self.coefRing.zero}
        if len(self.coefs)==0:
            self.coefs = {Monomial((0 for _ in range(len(self.vars))),self.vars): self.coefRing.zero}
    
    """Returns the multidegree of the polynomial."""
    def deg(self):
        monomials = self.coefs.keys()
        return (max([x.deg[i] for x in monomials]) for i in range(len(monomials[0].deg)))
    
    """Returns the leading term of the polynomial as a tuple."""
    def lt(self):
        k = next(iter(self.coefs))
        return (k,self.coefs[k])
    
    """Returns the leading monomial of the polynomial."""
    def lm(self):
        return next(iter(self.coefs))
    
    """Returns the leading coefficiente of the polynomial."""
    def lc(self):
        return self.coefs[next(iter(self.coefs))]
    
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
    
    """Product of two polynomials."""
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

    def __truediv__(self,other):
        return div_poly_RNF(self,other)[0][0]

    def __mod__(self, other):
        return div_poly_RNF(self,other)[1]

                
    def __eq__(self,other):
        # Monomials match and have same coefficients
        return hasattr(other, "val") and self.coefs == other.coefs
    
    def __str__(self):

        if self == type(self).zero:
            return "0"
        return " + ".join([f"{v}{'' if k.is_zero() else k}" for k,v in self.coefs.items() if v != type(v).zero])
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
        return hash((type(self).__name__, str(self.coefs)))

# Polynomials over integral domains
"""
Class representing a multivariate polynomial ring.
"""
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
        
    def __hash__(cls):
        return hash((type(cls).__name__, type(cls.coefRing).__name__, str(cls.vars)))
        
@external
def GetMultiPoly(ring, vars):
    attrs = {"coefRing": ring, "vars": vars}
    return MultivariatePolynomialRing(f"{ring}[{vars}]", (MultivariatePolynomial,), attrs)
