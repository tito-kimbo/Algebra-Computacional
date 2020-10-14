from rings import *

# Univariate polynomials over a Ring
class Polynomial:

    def __init__(self,coefs):
        self.coefs = coefs
    
    def deg(self):
        return len(self.coefs)-1
    
    # ziplongest would be more simple ut zero element is not accessible
    def __add__(self,other):
        d = min(self.deg(),other.deg())+1
        D = max(self.deg(),other.deg())+1
        aux = [self.coefs[i]+other.coefs[i] for i in range(d)]
        for i in range(d,self.deg()):
            aux.append(self.coefs[i])
        for i in range(d,other.deg()):
            aux.append(other.coefs[i])
        return Polynomial(aux)

    def __sub__(self,other):
        d = min(self.deg(),other.deg())+1
        D = max(self.deg(),other.deg())+1
        aux = [self.coefs[i]-other.coefs[i] for i in range(d)]
        for i in range(d,self.deg()):
            aux.append(self.coefs[i])
        for i in range(d,other.deg()):
            aux.append(other.coefs[i].opp())
        return Polynomial(aux)
    
    def __mul__(self,other):
        # TBD
        D = self.deg()*other.deg()+1
        aux = [None for _ in range(D)]
                
    def __eq__(self,other):
        return self.coefs == other.coefs
    
    def __str__(self):
        s = str(self.coefs[0])
        for i in range(1,self.deg()):
            s += " + " + str(self.coefs[i]) + " X^" + str(i)

# Polynomials over integral domains, WIP
class PolynomialRing(IntegralDomain):
        
    def __init__(self,ring,zero,one,elementClass=None):
        assert isinstance(ring,Ring)
        if elementClass is None:
            elementClass = Polynomial
        self.ring = ring
        super().__init__(zero,one,elementClass)