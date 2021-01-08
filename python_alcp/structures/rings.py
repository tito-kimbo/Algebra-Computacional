from abc import abstractmethod
from math import floor, sqrt
import random

from python_alcp.utils import (
        assuming,
        fast_exponentiation,
        repeated_addition,
        op_typecheck,
        externals
)

class Ring(type):

    def __init__(cls, name, bases, attrs):
        super().__init__(name, (*bases, RingElement), attrs)
        cls.zero = cls(0)
        cls.one = cls(1)
        cls.ring = cls
        cls.Element = cls

    def __eq__(cls, other):
        return cls.__name__ == other.__name__
    
    def __str__(cls):
        return cls.__name__
    
    def __repr__(cls):
        return Ring.__str__(cls)

    def __truediv__(cls, other):
        return externals.GetQuotient(cls, other)

    def __mul__(cls, other):
        return externals.GetIdeal(cls, [other])

    def __rmul__(cls, other):
        return Ring.__mul__(cls, other)

    def __getitem__(cls, key):
        return externals.GetPolynomials(cls, key)

    repr = "reduced"
    def setRepr(cls, rep):
        cls.repr = rep

    def is_integral(cls):
        return False

    def is_ufd(cls):
        return False

    def is_euclidean(cls):
        return False

    def is_field(cls):
        return False

    def is_polynomial(cls):
        return False

    @abstractmethod
    def char(cls):
        pass

    @abstractmethod
    def order(cls):
        pass

    @abstractmethod
    def is_finite(cls):
        pass

    @abstractmethod
    def generators(cls):
        pass

    def elements(cls):
        # Note: potentially infinite!

        if not hasattr(cls, "__exploring"):
            cls.__exploring = {cls.one}.union(cls.generators())
            cls.__visited = {cls.zero}

        def elemgen(exploring, visited):

            for e in visited:
                yield e

            while len(exploring) > 0:
                e1 = random.sample(exploring, 1)[0]
                if e1 not in visited:
                    yield e1
                visited.add(e1)
                if -e1 not in visited:
                    exploring.add(-e1)
                for e2 in visited:
                    if (e1+e2) not in visited:
                        exploring.add(e1+e2)
                    if (e1*e2) not in visited:
                        exploring.add(e1*e2)
                exploring.remove(e1)

        return elemgen(cls.__exploring, cls.__visited)

    @abstractmethod
    def units(cls):
        pass

    def build(cls, *args, **kwargs):
        return cls(*args, **kwargs)

    def __hash__(cls):
        return hash(cls.__name__)


class IntegralDomain(Ring):

    def __init__(cls, name, bases, attrs):
        return super().__init__(name, (*bases, IntegralDomainElement), attrs)

    def is_integral(cls):
        return True

class UniqueFactorizationDomain(IntegralDomain):

    def __init__(cls, name, bases, attrs):
        return super().__init__(name, bases, attrs)

    def is_ufd(cls):
        return True


class EuclideanDomain(UniqueFactorizationDomain):

    def __init__(cls, name, bases, attrs):
        return super().__init__(name, (*bases, EuclideanDomainElement), attrs)

    def is_euclidean(cls):
        return True

    def phi(self,element):
        pass

    def numphi(self, n):
        pass

class Field(EuclideanDomain):

    def phi(self,element):
        return 1;

    def numphi(self, n):
        return 0 if n <= 1 else self.char()

    def is_field(cls):
        return True


class RingElement():

    val = None

    def __init__(self, val):
        self.val = val

    def __add__(self,other):
        pass

    def __sub__(self,other):
        pass

    def inner_mul(self, other):
        """ Multiplication of ring elements """
        pass

    def __mul__(self,other):
        if isinstance(other, int):
            return repeated_addition(self, other)
        elif isinstance(other,Ring):
            return other.__mul__(self)
        else:
            return self.inner_mul(other)
        pass

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        assuming(other.is_unit(), "Can't divide by non-unit")
        return self * other.inverse()

    def inverse(self):
        assuming(self.is_unit(), "Can't invert non-unit")
        return next(u for u in type(self).units() if u*self == type(self).one)

    def __pow__(self, other):
        if type(other) == int:
            if other < 0:
                return fast_exponentiation(self.inverse(), -other)
            return fast_exponentiation(self, other)

    @abstractmethod
    def is_unit(self):
        pass

    @abstractmethod
    def is_prime(self):
        pass

    def __eq__(self,other):
        return hasattr(other, "val") and self.val == other.val

    @abstractmethod
    def __str__(self):
        pass
    
    def __repr__(self):
        return self.__str__()

    @abstractmethod
    def __neg__(self):
        pass

    def __hash__(self):
        return hash((type(self).__name__, self.val))

    def __lt__(self,other):
        pass

    def __le__(self, other):
        return self.__lt__(other) or self == other

    def __gt__(self, other):
        return not self.__le__(other)

    def __ge__(self, other):
        return not self.__lt__(other)

    def normal(self):
        if self.is_unit():
            return type(self).one
        else:
            return self

    def __int__(self):
        val = self.val
        while hasattr(val, "val"):
            val = val.val
        return int(val)


class IntegralDomainElement(RingElement):

    @abstractmethod
    def __floordiv__(self,other):
        pass

    @abstractmethod
    def __mod__(self,other):
        pass

    def in_base(self, b):
        """
            returns the element expressed in base b
        """

        R = self.ring
        x = self
        res = []

        while x != R.zero:
            res.append(x % b)
            x = x // b

        return res[::-1]



class UniqueFactorizationDomainElement(IntegralDomainElement):

    @abstractmethod
    def factors(self):
        pass


class EuclideanDomainElement(UniqueFactorizationDomainElement):
    pass



class FieldElement(EuclideanDomainElement):

    @abstractmethod
    def inverse(self):
        pass

    @abstractmethod
    def __truediv__(self,other):
        pass

    def is_unit(self):
        return True

    def is_prime(self):
        return True

    def factors(self):
        return [self]

    def normal(self):
        return type(self).one
