from math import floor, sqrt

from python_alcp.utils import assuming, try_op, external, primes, prime_factors
from python_alcp.structures.rings import EuclideanDomain, EuclideanDomainElement


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
        return hash(self.val)

    def __lt__(self,other):
        while hasattr(other, 'val'):
            other = other.val
        return self.val < other

    def is_prime(self):
        for p in primes():
            if self.val == p:
                return True
            if p > self:
                return False

    def is_unit(self):
        return self.val == 1 or self.val == -1

    def factors(self):
        return {type(self)(k):v for k,v in prime_factors(self).items()}




class Integers(EuclideanDomain):

    def units(cls):
        return {cls(-1), cls(1)}

    def generators(cls):
        return {cls(1)}

    def phi(cls, element):
        return abs(element.val)

    def numphi(cls, n):
        return n



Z = Integers("Z", (IntElement,), {})
external(Z)
