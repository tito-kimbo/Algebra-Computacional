from structures.rings import *
from itertools import zip_longest
from utils import assuming

VARS = ["X","Y","Z","T","U","V"] # Could be extended arbitrarily with sub indexing

def addAll(l,neutral):
    res = neutral
    for x in l:
        res = res+x
    return res

def multAll(l,neutral):
    res = neutral
    for x in l:
        res = res*x
    return res

# Polynomials over integral domains, WIP
class PolynomialRing(IntegralDomain):
        
    def __init__(self,ring,elementClass=None):
        assuming(isinstance(ring,Ring), "ring must be a Ring")
            
        # For string conversion
        self.chain = 0
        if isinstance(ring,PolynomialRing):
            self.chain = ring.chain+1
        
        chain = self.chain
        # Univariate polynomials over a Ring - dynamically linked to the current ring
        class Polynomial:

            def __init__(self,coefs):
                for c in coefs:
                    assuming(type(c)==ring.elementClass, "coefficients must be ring elements")
                self.coefs = coefs
            
            def deg(self):
                return len(self.coefs)-1
            
            def __add__(self,other):
                return Polynomial([x+y for x,y in zip_longest(self.coefs,other.coefs,fillvalue=ring.zero)])

            def __sub__(self,other):
                return Polynomial([x-y for x,y in zip_longest(self.coefs,other.coefs,fillvalue=ring.zero)])
            
            # Convolutional product
            def __mul__(self,other):
                D = self.deg()+other.deg()+1
                c,o = self.coefs[:],other.coefs[:] # shallow copies to avoid issues with original
                while len(c) < D:
                    c.append(ring.zero)
                coefs = list(zip_longest(c,o,fillvalue=ring.zero))
                aux = []
                for i in range(D):
                    aux.append(addAll([coefs[k][0]*coefs[i-k][1] for k in range(i+1)],ring.zero))
                while len(aux)!=0 and aux[-1] == ring.zero:
                    aux.pop(-1)
                return Polynomial(aux)
                        
            def __eq__(self,other):
                return self.coefs == other.coefs
            
            def __str__(self):
                s = "(" + str(self.coefs[0])
                for i in range(1,self.deg()+1):
                    s += " + " + str(self.coefs[i]) + "*" + VARS[chain] + "^" + str(i)
                s += ")"
                return s
        
        if elementClass is None:
            elementClass = Polynomial
        self.ring = ring
        super().__init__(Polynomial([ring.zero]),Polynomial([ring.one]),elementClass)
    
    def build(self,args):
        l = []
        for i in range(len(args)):
            l.append(self.ring.build(args[i]))
        return self.elementClass(l)
        
