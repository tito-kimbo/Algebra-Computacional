from structures.rings import EuclideanDomain
from math import floor, sqrt
from utils import assuming


OPERAND_ERROR = "All operands must be within the following list "
def _op_typecheck(operand,allowed):
    if not type(operand) in allowed:
        raise TypeError(OPERAND_ERROR+" "+str(allowed))


"""
    Ring of gaussian integers
    Usage:
        from structures.rings import Zi
        gi_i = Zi.build(0,1)
        gi2 = Zi.one + Zi.one
"""
class Z(EuclideanDomain):
    
    class Element(EuclideanDomain.Element):
    
        def __init__(self,val, *args, **kw):
            super().__init__(*args, **kw)
            _op_typecheck(val,allowed=[int])
            self.val=val
        
        def __add__(self,other):
            _op_typecheck(other,allowed=[Z.Element])
            return Z.Element(self.val+other.val)

        def __sub__(self,other):
            _op_typecheck(other,allowed=[Z.Element])
            return Z.Element(self.val-other.val)
        
        def inner_mul(self,other):
            _op_typecheck(other,allowed=[Z.Element])
            return Z.Element(self.val*other.val)
        
        def __eq__(self,other):
            return type(other) is type(self) and self.val==other.val
        
        def __floordiv__(self,other):
            _op_typecheck(other,allowed=[Z.Element])
            return Z.Element(self.val//other.val)

        def __truediv__(self,other):
            _op_typecheck(other,allowed=[Z.Element])
            return Z.Element(self.val//other.val)

        def __mod__(self,other):
            _op_typecheck(other,allowed=[Z.Element])
            return Z.Element(self.val%other.val)
            
        def __neg__(self):
            return Z.Element(-self.val)

        def normal(self):
            if(self.val < 0):
                return -self
            else:
                return self 

        def __str__(self):
            return str(self.val)

        def __hash__(self):
            return hash(("Z", self.val))

        def __lt__(self,other):
            _op_typecheck(other,allowed=[Z.Element])
            return self.val < other.val

        def is_prime(self):
            # Temporal
            for i in range(2,floor(sqrt(self.val)) + 1):
                if self.val % i == 0:
                    return False
            # TODO negative values??
            return self.val > 1

        def is_unit(self):
            return self.val == 1 or self.val == -1

        def factors(self):

            def add(d,f):
                if f in d:
                    d[f] += 1
                else:
                    d[f] = 1

            # Temporal
            val = self.val
            res = dict()
            while(val > 1):
                found = False
                for i in range(2,floor(sqrt(val)) + 1):
                    if val % i == 0:
                        found = True
                        add(res,i)
                        val = val//i
                        break
                if not found:
                    add(res,val)
                    break
            # TODO negative values??
            return res

            
    def __init__(self):
        super().__init__(self.build(0),self.build(1))
    
    def phi(self,element):
        if type(element) is not self.Element:
            raise TypeError("Phi can only be applied to elements of the ring")
        return abs(element.val)

    def numphi(self, n):
        return n

    def __eq__(self, other):
        return other.__class__ == self.__class__

    def __str__(self):
        return "\N{DOUBLE-STRUCK CAPITAL Z}"

    
    def __mul__(self, other):
        if type(other) == int:
            return super().__mul__(self.build(other))
        else:
            return super().__mul__(other)


    def char(self):
        return 0

    def order(self):
        return -1


Z = Z()

"""
    Ring of integers
    Usage:
        from structures.rings import Z
        i18 = Z.build(18)
        i2 = Z.one + Z.one
"""
class GaussianIntegers(EuclideanDomain):
    
    class Element(EuclideanDomain.Element):
    
        def __init__(self,a,b, *args, **kw):
            super().__init__(*args, **kw)
            _op_typecheck(a,allowed=[int])
            _op_typecheck(b,allowed=[int])
            self.a=a
            self.b=b
        
        def __add__(self,other):
            _op_typecheck(other,allowed=[Zi.Element])
            return Zi.Element(self.a+other.a,self.b+other.b)

        def __sub__(self,other):
            _op_typecheck(other,allowed=[Zi.Element])
            return Zi.Element(self.a-other.a,self.b-other.b)
        
        def inner_mul(self,other):
            _op_typecheck(other,allowed=[Zi.Element])
            return Zi.Element(self.a*other.a-self.b*other.b,self.a*other.b+self.b*other.a)
        
        def __eq__(self,other):
            return type(other) is type(self) and self.a==other.a and self.b==other.b
        
        # WIP
        def __floordiv__(self,other):
            _op_typecheck(other,allowed=[Zi.Element])
            return Zi.Element(self.val//other.val)
        # WIP
        def __truediv__(self,other):
            _op_typecheck(other,allowed=[Zi.Element])
            return Zi.Element(self.val//other.val)
        # WIP
        def __mod__(self,other):
            _op_typecheck(other,allowed=[Zi.Element])
            return Zi.Element(self.val%other.val)
            
        def __neg__(self):
            return Zi.Element(-self.a,-self.b)
        
        def conj(self):
            return Zi.Element(self.a,-self.b)

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
            return hash(("Z[i]", self.a, self.b))
        
        def is_unit(self):
            return (self.a == 0 and self.b in [-1,1]) or (self.a in [-1,1] and self.b == 0)
            
        def factors(self,other):
            pass
    
        def is_prime(self,other):
            pass
        
        def __lt__(self,other):
            raise NotImplementedError("This should never be called.")


            
    def __init__(self):
        super().__init__(self.build(0,0),self.build(1,0))
    
    def phi(self,element):
        if type(element) is not self.Element:
            raise TypeError("Phi can only be applied to elements of the ring")
        return sqrt(element.a*element.a+element.b*element.b)

    def numphi(self, n):
        return n

    def __eq__(self, other):
        return other.__class__ == self.__class__

    def __str__(self):
        return "\N{DOUBLE-STRUCK CAPITAL Z}[i]"

    
    def __mul__(self, other):
        # Assumes other is of type GaussianIntegers.Element
        if type(other) == GaussianIntegers.Element:
            return super().__mul__(other)
        else:
            raise ValueError("Incompatible type")
    
    def char(self):
        return 0

    def order(self):
        return -1

Zi = GaussianIntegers()