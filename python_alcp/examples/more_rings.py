from math import floor, sqrt

from python_alcp.utils import (
    assuming,
    try_op,
    external,
    int_gcd,
)
from python_alcp.structures.rings import Field, FieldElement



######### Rationals ###############


class RationalElement(FieldElement):

    def __init__(self,num,den = 1):

        while hasattr(num, "val"):
            num = num.val

        if hasattr(num, "__len__") and len(num) == 2:
            num, den = num

        assuming(hasattr(num, "__int__") and hasattr(den, "__int__") and int(den) != 0, 
                f"Can't build a rational from {num}, {den}")

        num, den = int(num), int(den)

        g = int_gcd(num,den)

        self.num = num//g
        self.den = den//g
        self.val = (self.num, self.den)

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return f"{self.num}/{self.den}"

    @try_op()
    def __add__(self,other):
        if not hasattr(other, "den"):
            other = type(self)(other)
        return type(self)(self.num*other.den + self.den*other.num, self.den*other.den)

    @try_op()
    def __sub__(self,other):
        if not hasattr(other, "den"):
            other = type(self)(other)
        return type(self)(self.num*other.den - self.den*other.num, self.den*other.den)
    
    @try_op("multiply")
    def inner_mul(self,other):
        if not hasattr(other, "den"):
            other = type(self)(other)
        return type(self)(self.num*other.num, self.den*other.den)
    
    def __eq__(self,other):
        return isinstance(other, type(self)) and self.val==other.val

    @try_op("divide")
    def __truediv__(self,other):
        if not hasattr(other, "den"):
            other = type(self)(other)
        return type(self)(self.num*other.den, self.den*other.num)
        
    def __neg__(self):
        return type(self)(-self.num, self.den)

    def normal(self):
        if(self.val < 0):
            return -self
        else:
            return self 

    def __hash__(self):
        return hash((type(self).__name__, self.val))

    def __lt__(self,other):
        if not hasattr(other, 'den'):
            other = type(self)(other)
        return self.num / self.den < other.num / other.den

    def inverse(self):
        return type(self)(self.den, self.num)




class Rationals(Field):

    def units(cls):
        return {}

    def generators(cls):
        return {cls(1)}

    def elements(cls):
        def elemgen():
            yield 0
            i = 1
            while True:
                j = 1
                while j < i:
                    yield cls(j,i-j)
                    yield cls(-j,i-j)
                    j += 1;
                i += 1
        return elemgen()

    def phi(cls, element):
        return 1

    def numphi(cls, n):
        return 0


Q = Rationals("\N{DOUBLE-STRUCK CAPITAL Q}", (RationalElement,), {})
external(Q, name="Q")
