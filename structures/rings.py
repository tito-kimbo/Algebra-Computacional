from abc import ABC,abstractmethod
from utils import assuming

class Ring(ABC):
    """Class representing an algebraic ring (commutative ring with unity)"""
    
    class Element(ABC): 
        """Class representing an element of the ring."""

        ring = None

        @abstractmethod
        def __add__(self,other):
            assuming(self.ring == other.ring,
                    "You can only add elements of the same ring")
            pass

        @abstractmethod
        def __sub__(self,other):
            assuming(self.ring == other.ring,
                    "You can only substract elements of the same ring")
            pass

        @abstractmethod
        def inner_mul(self, other):
            """ Multiplication of ring elements """
            pass

        def __mul__(self,other):
            if type(other) == int:
                return repeated_addition(self, other)
            elif isinstance(other,Ring):
                return other.__mul__(self)
            else:
                return self.inner_mul(other)
            pass

        def __rmul__(self, other):
            return self.__mul__(other)

        def __pow__(self, other):
            if type(other) == int:
                return fast_exponentiation(self, other)

        @abstractmethod
        def is_unit(self):
            pass

        @abstractmethod
        def __eq__(self,other):
            pass
            
        @abstractmethod
        def __str__(self):
            pass
        
        def __repr__(self):
            return self.__str__()

        @abstractmethod
        def __neg__(self):
            pass

        @abstractmethod
        def __hash__(self):
            pass

        def __le__(self, other):
            return self.__lt__(other) or self == other

        def __gt__(self, other):
            return not self.__le__(other)

        def __ge__(self, other):
            return not self.__lt__(other)


    def __new__(cls, *args, **kwargs): 
        # Decorates the Element class before construction of the ring
        # so that self.Element.ring == self

        obj = super(Ring, cls).__new__(cls)

        obj.Element = type('.'.join([cls.__module__,cls.Element.__qualname__]), (cls.Element,), {"ring": obj})
        return obj

    def __init__(self,zero,one, **kw):
        self.zero = zero
        self.one = one

        super().__init__(**kw)
    
    def build(self,*args,**kwargs):
        return self.Element(*args, **kwargs)

    @abstractmethod
    def __eq__(self, other):
        pass
    
    @abstractmethod
    def __str__(self):
        pass
    
    def __repr__(self):
        return self.__str__()

    # Adds support for stuff like Z/NZ(5) instead of Quotient(Z, NZ(5))
    def __truediv__(self, other):
        from structures.ideals import GetQuotient
        return GetQuotient(self, other)

    # Adds support for 6*Z
    def __mul__(self, other):
        assuming(other.ring == self, f"Can't generate ideal from {other} in {self}")
        from structures.ideals import GetIdeal
        return GetIdeal([other])

    def __rmul__(self, other):
        return self.__mul__(other)

    # Adds support for Z['x']
    def __getitem__(self, key):
        from structures.polynomials import GetPolynomials
        return GetPolynomials(self, key)
    
    def is_integral(self):
        return False

    def is_euclidean(self):
        return False

    @abstractmethod
    def char(self):
        pass

    @abstractmethod
    def order(self):
        pass

class IntegralDomain(Ring):
    """Class representing an integral domain."""

    class Element(Ring.Element): 
        """Class representing an element of the ID."""
        @abstractmethod
        def __floordiv__(self,other):
            assuming(self.ring == other.ring,
                    "You can only divide elements of the same ring")
            pass
        @abstractmethod
        def __mod__(self,other):
            assuming(self.ring == other.ring,
                    "You can only calculate the remainder between elements of the same ring")
            pass

        @abstractmethod
        def is_prime(self):
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

    def is_integral(self):
        return True


class UniqueFactorizationDomain(IntegralDomain):

    class Element(IntegralDomain.Element):

        @abstractmethod
        def factors(self):
            pass

    
class EuclideanDomain(UniqueFactorizationDomain):
    """Structure representing an euclidean domain."""

    def is_euclidean(self):
        return True
    
    @abstractmethod
    def phi(self,element):
        pass

    @abstractmethod
    def numphi(self, n: int):
        """ 
            Returns the number of elements e such that phi(e) < n, or -1
            if the number is infinite / unknown
        """
        pass


def repeated_addition(element: Ring.Element, n: int):
    """ Computes a+a+a+a... (n times) in O(log(n)) time """

    res = element.ring.zero
    acc = element
    while n > 0:
        if n % 2 == 1:
            res += acc
        n //= 2
        acc = acc+acc

    return res

def fast_exponentiation(element: Ring.Element, n: int):
    """
        Computes a^n in O(log(n)) time
    """

    res = element.ring.one
    acc = element
    while n > 0:
        if n % 2 == 1:
            res *= acc
        n //= 2
        acc = acc*acc

    return res
