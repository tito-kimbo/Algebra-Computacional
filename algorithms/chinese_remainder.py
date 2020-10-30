from structures.ideals import *
from algorithms.divisibility import gcd
from functools import reduce

# WIP
def modinv(x,f,id):
    """Receives two elements of the same ring R and returns the inverse of x in R/<f>. <f> must be maximal.
       
       :x  element of R
       :f  element of R generating a maximal ideal
       :id ideal class """
    if x.ring != f.ring:
        raise ValueError('Both elements must belong to the same ring')
    R = x.ring
    Q = R // id(f)
    return Q.build(x.val).inverse().val 
    

def chinese_remainder(eqs,mods,id):
    """Receives as argument a list of modular equations and returns the solution. The moduli must be pairwise coprime.
       
       :eqs    the ith element represents the value of x modulo the ith element of mods
       :mods   the modulo for each equation
       :id     ideal class"""
    if len(eqs) != len(mods):
        raise Exception('Not a proper equation representation')
    L = len(eqs)
    for i,j in zip(range(L),range(L)):
        if i != j and gcd(mods[i].mods[j]) != mods[i].ring.one():
            raise ValueError('The moduli are not pairwise coprime')

    N = reduce((lambda x,y:x*y), mods)
    r = reduce((lambda x,y:x+y), [eqs[i]*(N/mods[i])*modinv(N/mods[i],mods[i],id) for i in range(L)])
    return r
