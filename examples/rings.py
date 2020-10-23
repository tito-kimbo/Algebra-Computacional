from structures.rings import *

OPERAND_ERROR = "All operands must be within the following list "
def _op_typecheck(operand,allowed):
            if not type(operand) in allowed:
                raise TypeError(OPERAND_ERROR+" "+str(allowed))


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



class ZModFactory:
    """Utility class that produces any ring of integers modulo n"""

    # Element methods
        
    def _elem_init(self,val):
        _op_typecheck(val,allowed=[int])
        self.val = val % self.mod
        
    def _elem_add(self,other):
        _op_typecheck(other,allowed=[self.__class__])
        return self.__class__(self.val+other.val)

    def _elem_sub(self,other):
        _op_typecheck(other,allowed=[self.__class__])
        return self.__class__(self.val-other.val)
    
    def _elem_mul(self,other):
        _op_typecheck(other,allowed=[self.__class__])
        return self.__class__(self.val*other.val)
    
    def _elem_eq(self,other):
        return type(other) is type(self) and self.val==other.val
    
    def _elem_truediv(self,other):
        _op_typecheck(other,allowed=[self.__class__])
        return self.__class__(self.val//other.val)
    
    def _elem_mod(self,other):
        _op_typecheck(other,allowed=[self.__class__])
        return self.__class__(self.val%other.val)
    
    def _elem_str(self):
        return str(self.val)

    # Structure methods

    def _mod_init(self, elementClass):
        super(EuclideanDomain, self).__init__(elementClass(0),elementClass(1),elementClass)

    def _phi(self,element):
        if type(element) is not self.elementClass:
            raise TypeError("Phi can only be applied to elements of the ring")
        return abs(element.val)



    def create_zmod(mod: int):
        """ Returns a class representing the integers modulo mod
            Example usage:
                Z5 = ZModFactory.create_zmod(5)
                struct = Z5()
                i2 = Z5.SymbolicModularInteger(2)
                i3 = Z5.SymbolicModularInteger(3)
                (i2+i3) == struct.zero                # True
                (struct.zero - struct.one).val == 5   # True
        """

        if mod < 2:
            raise ValueError(f"Can't create Zmod, mod must be >= 2, received {mod}")

        elementClass = type(f"Z{mod}.SymbolicModularInteger", (EuclideanDomain.IDElement, ),
                {"mod": mod,
                 "__init__": ZModFactory._elem_init,
                 "__add__": ZModFactory._elem_add,
                 "__sub__": ZModFactory._elem_sub,
                 "__mul__": ZModFactory._elem_mul,
                 "__eq__": ZModFactory._elem_eq,
                 "__truediv__": ZModFactory._elem_truediv,
                 "__mod__": ZModFactory._elem_mod,
                 "__str__": ZModFactory._elem_str})

        fieldClass = type(f"Z{mod}", (EuclideanDomain, ),
                {"mod": mod,
                 "__init__": lambda x: ZModFactory._mod_init(x, elementClass),
                 "SymbolicModularInteger": elementClass,
                 "elementClass": elementClass,
                 "phi": ZModFactory._phi})

        return fieldClass
