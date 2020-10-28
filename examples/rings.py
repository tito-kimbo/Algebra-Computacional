from structures.rings import EuclideanDomain
from structures.ideals import Ideal
from math import floor, sqrt


OPERAND_ERROR = "All operands must be within the following list "
def _op_typecheck(operand,allowed):
    if not type(operand) in allowed:
        raise TypeError(OPERAND_ERROR+" "+str(allowed))


"""
    Ring of integers
    Usage:
        from structures.rings import Z
        i18 = Z.build(18)
        i2 = Z.one + Z.one
"""
class Z(EuclideanDomain):
    
    class Element(EuclideanDomain.Element):
    
        def __init__(self,val):
            _op_typecheck(val,allowed=[int])
            self.val=val
        
        def __add__(self,other):
            _op_typecheck(other,allowed=[Z.Element])
            return Z.Element(self.val+other.val)

        def __sub__(self,other):
            _op_typecheck(other,allowed=[Z.Element])
            return Z.Element(self.val-other.val)
        
        def __mul__(self,other):
            _op_typecheck(other,allowed=[Z.Element])
            return Z.Element(self.val*other.val)
        
        def __eq__(self,other):
            return type(other) is type(self) and self.val==other.val
        
        def __truediv__(self,other):
            _op_typecheck(other,allowed=[Z.Element])
            return Z.Element(self.val//other.val)
        
        def __mod__(self,other):
            _op_typecheck(other,allowed=[Z.Element])
            return Z.Element(self.val%other.val)
            
        def opp(self,other):
            _op_typecheck(other,allowed=[Z.Element])
            return Z.symbolicInteger(-1)*other

        def __str__(self):
            return str(self.val)

        def _is_prime(self):
            # Temporal
            for i in range(2,floor(sqrt(self.val)) + 1):
                if self.val % i == 0:
                    return False
            return True
            
    def __init__(self):
        super().__init__(Z.Element(0),Z.Element(1))
    
    def phi(self,element):
        if type(element) is not Z.Element:
            raise TypeError("Phi can only be applied to elements of the ring")
        return abs(element.val)


class NZ(Ideal):
    def __init__(self,generator):
        if isinstance(generator, int):
            generator = Z.build(generator)
        super().__init__(Z,[generator])

        # Is there a better way? A general one?
        self._maximal = generator._is_prime()
    
    def has(self,element):
        return element % self.generators[0] == self.ring.zero

    def is_principal(self):
        return True

    def is_maximal(self):
        return self._maximal


Z = Z()
