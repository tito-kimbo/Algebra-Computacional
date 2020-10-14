from rings import *
from itertools import zip_longest

VARS = ["X","Y","Z","T","U","V"] # Could be extended arbitrarily with sub indexing

# Polynomials over integral domains, WIP
class PolynomialRing(IntegralDomain):
        
    def __init__(self,ring,zero,one,elementClass=None):
        assert isinstance(ring,Ring)
        
        # For string conversion
        self.chain = 0
        if isinstance(ring,PolynomialRing):
            self.chain = ring.chain+1
        
        # Univariate polynomials over a Ring - dynamically linked to the current ring
        class Polynomial:

            def __init__(self,coefs):
                self.coefs = coefs
            
            def deg(self):
                return len(self.coefs)-1
            
            def __add__(self,other):
                return Polynomial([x+y for x,y in zip_longest(self.coefs,other.coefs,fillvalue=zero)])

            def __sub__(self,other):
                return Polynomial([x-y for x,y in zip_longest(self.coefs,other.coefs,fillvalue=zero)])
            
            # Convolutional product
            def __mul__(self,other):
                D = self.deg()*other.deg()+1
                coefs = zip_longest(self.coefs,other.coefs,fillvalue=one)
                aux = []
                for i in range(D):
                    aux.append(sum([coefs[k][0]*coefs[i-k][1] for k in range(i)]))
                return Polynomial(aux)
                        
            def __eq__(self,other):
                return self.coefs == other.coefs
            
            def __str__(self):
                s = "(" + str(self.coefs[0])
                for i in range(1,self.deg()):
                    s += " + " + str(self.coefs[i]) + "*" + VARS[self.chain] + "^" + str(i)
                s += ")"
            
        if elementClass is None:
            elementClass = Polynomial
        self.ring = ring
        super().__init__(zero,one,elementClass)