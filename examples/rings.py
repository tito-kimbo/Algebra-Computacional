from structures.rings import EuclideanDomain
from structures.ideals import Ideal


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
    
    class SymbolicInteger(EuclideanDomain.IDElement):
    
        def __init__(self,val):
            _op_typecheck(val,allowed=[int])
            self.val=val
        
        def __add__(self,other):
            _op_typecheck(other,allowed=[Z.SymbolicInteger])
            return Z.SymbolicInteger(self.val+other.val)

        def __sub__(self,other):
            _op_typecheck(other,allowed=[Z.SymbolicInteger])
            return Z.SymbolicInteger(self.val-other.val)
        
        def __mul__(self,other):
            _op_typecheck(other,allowed=[Z.SymbolicInteger])
            return Z.SymbolicInteger(self.val*other.val)
        
        def __eq__(self,other):
            return type(other) is type(self) and self.val==other.val
        
        def __truediv__(self,other):
            _op_typecheck(other,allowed=[Z.SymbolicInteger])
            return Z.SymbolicInteger(self.val//other.val)
        
        def __mod__(self,other):
            _op_typecheck(other,allowed=[Z.SymbolicInteger])
            return Z.SymbolicInteger(self.val%other.val)
            
        def opp(self,other):
            _op_typecheck(other,allowed=[Z.SymbolicInteger])
            return Z.symbolicInteger(-1)*other
        
        def __str__(self):
            return str(self.val)
            
    def __init__(self):
        super().__init__(Z.SymbolicInteger(0),Z.SymbolicInteger(1),Z.SymbolicInteger)
    
    def phi(self,element):
        if type(element) is not Z.SymbolicInteger:
            raise TypeError("Phi can only be applied to elements of the ring")
        return abs(element.val)

class NZ(Ideal):
    def __init__(self,generator):
        super().__init__(Z(),[generator])
    
    def has(self,element):
        return element % self.generators[0] == self.ring.zero

Z = Z()
