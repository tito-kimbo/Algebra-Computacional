from math import floor, sqrt

from python_alcp.utils import (
    assuming,
    try_op,
    external,
    primes,
    prime_factors,
    op_typecheck
)
from python_alcp.structures.rings import EuclideanDomain, EuclideanDomainElement



######### Integers ###############

"""
    Ring of integers
    Usage:
        from python_alcp.examples.rings import Z
        i18 = Z(18)
        i2 = Z.one + Z.one
"""

class IntElement(EuclideanDomainElement):

    def __init__(self,val):

        while hasattr(val, "val"):
            val = val.val

        assuming(hasattr(val, "__int__"), f"Can't build an integer from {val}")
        self.val = int(val)

    def __repr__(self):
        return str(self.val)

    def __str__(self):
        return str(self.val)

    @try_op()
    def __add__(self,other):
        while hasattr(other, 'val'):
            other = other.val
        return type(self)(self.val + other)

    @try_op()
    def __sub__(self,other):
        while hasattr(other, 'val'):
            other = other.val
        return type(self)(self.val - other)
    
    @try_op("multiply")
    def inner_mul(self,other):
        while hasattr(other, 'val'):
            other = other.val
        return type(self)(self.val * other)
    
    def __eq__(self,other):
        return isinstance(other, type(self)) and self.val==other.val

    @try_op("divide")
    def __floordiv__(self,other):
        while hasattr(other, 'val'):
            other = other.val
        return type(self)(self.val // other)

    @try_op("divide")
    def __truediv__(self,other):
        while hasattr(other, 'val'):
            other = other.val
        return type(self)(self.val / other)

    @try_op("calculate remainder")
    def __mod__(self,other):
        while hasattr(other, 'val'):
            other = other.val
        return type(self)(self.val % other)
        
    def __neg__(self):
        return type(self)(-self.val)

    def normal(self):
        if(self.val < 0):
            return -self
        else:
            return self 

    def __str__(self):
        return str(self.val)

    def __hash__(self):
        return hash((type(self).__name__, self.val))

    def __lt__(self,other):
        while hasattr(other, 'val'):
            other = other.val
        return self.val < other

    def is_prime(self):
        for p in primes():
            if self.val % p == 0:
                return False
            if p**2 > self:
                return True

    def is_unit(self):
        return self.val == 1 or self.val == -1

    def inverse(self):
        assuming(self.is_unit(), "Can't invert non-unit")
        return self

    def factors(self):
        return {type(self)(k):v for k,v in prime_factors(self).items()}




class Integers(EuclideanDomain):

    def units(cls):
        return {cls(-1), cls(1)}

    def generators(cls):
        return {cls(1)}

    def elements(cls):
        def elemgen():
            yield 0
            i = 1
            while True:
                yield cls(i)
                yield cls(-i)
                i += 1
        return elemgen()

    def phi(cls, element):
        return abs(element.val)

    def numphi(cls, n):
        return n



Z = Integers("\N{DOUBLE-STRUCK CAPITAL Z}", (IntElement,), {})
external(Z, name="Z")






######### Gaussian Integers ###############



class GaussianIntegers(EuclideanDomain):
    
    def phi(cls,element):
        if not isinstance(element, cls):
            raise TypeError("Phi can only be applied to elements of the ring")
        return sqrt(element.a*element.a+element.b*element.b)

    def numphi(cls, n):
        return n

    def __mul__(cls, other):
        # Assumes other is of type GaussianIntegers
        if isinstance(other, cls):
            return super().__mul__(other)
        else:
            raise ValueError("Incompatible type")
    
    def char(cls):
        return 0

    def order(cls):
        return 0


class GaussianIntegerElement(EuclideanDomainElement):

    def __init__(self,a,b = None, *args, **kw):
        op_typecheck(a,allowed=[int])
        self.a=a
        if b is not None:
            op_typecheck(b,allowed=[int])
            self.b=b
        else:
            self.b = 0
    
    def __add__(self,other):
        op_typecheck(other,allowed=[type(self)])
        return type(self)(self.a+other.a,self.b+other.b)

    def __sub__(self,other):
        op_typecheck(other,allowed=[type(self)])
        return type(self)(self.a-other.a,self.b-other.b)
    
    def inner_mul(self,other):
        op_typecheck(other,allowed=[type(self)])
        return type(self)(self.a*other.a-self.b*other.b,self.a*other.b+self.b*other.a)
    
    def __eq__(self,other):
        return type(other) is type(self) and self.a==other.a and self.b==other.b
    
    def __floordiv__(self,other):
        op_typecheck(other,allowed=[type(self)])
        # Calculate quotient in the complex plane (via conjugate)
        frac = (self*other.conj(),other*other.conj())
        # Find sufficiently close element of Z[i]
        a1,a2 = frac[0].a//frac[1].a,frac[0].a//frac[1].a+1
        b1,b2 = frac[0].b//frac[1].a,frac[0].b//frac[1].a+1
        threshold = Zi.phi(other)
        
        for a in [a1,a2]:
            for b in [b1,b2]:
                c = Zi.build(a,b)
                if Zi.phi(self-other*c) < threshold:
                    return c
        raise RuntimeError("This should never happen.")
    
    def __truediv__(self,other):
        return (self//other)
    
    def __mod__(self,other):
        op_typecheck(other,allowed=[type(self)])
        return (self - other*(self/other))
        
    def __neg__(self):
        return type(self)(-self.a,-self.b)
    
    def conj(self):
        return type(self)(self.a,-self.b)

    def __str__(self):
        if self.b == 0:
            s = str(self.a)
        elif self.a == 0:
            s = ''.join([str(self.b),'i'])
        elif self.b < 0:
            s = ''.join([str(self.a),str(self.b),'i'])
        else:
            s = ''.join([str(self.a),'+',str(self.b),'i'])
        return ''.join(['(',s,')'])
    
    def __hash__(self):
        return hash((type(self).__name__, self.a, self.b))
    
    def is_unit(self):
        return (self.a == 0 and self.b in [-1,1]) or (self.a in [-1,1] and self.b == 0)
        
    def factors(self,other):
        pass

    def is_prime(self,other):
        pass
    
    def __lt__(self,other):
        raise NotImplementedError("This should never be called.")


Zi = GaussianIntegers("\N{DOUBLE-STRUCK CAPITAL Z}[i]", (GaussianIntegerElement,), {})
external(Zi, name="Zi")
