from structures.ideals import *
from algorithms.divisibility import gcd

# WIP
def modinv(x,f):
    """Receives two elements of the same ring R and returns the inverse of x in R/<f>. <f> must be maximal.
       
       :x Element of R
       :f Element of R generating a maximal ideal"""
    if x.ring != f.ring:
        raise ValueError('Both elements must belong to the same ring')
    R = x.ring
    Q = R // IDEAL(f) # Need to check how to implement this properly 
    return R.build(Q.build(x.value).inverse().value) # Requires testing
    

def chinese_remainder(eqs,mods):
    """Receives as argument a list of modular equations and returns the solution. The moduli must be pairwise coprime.
       
       :eqs    the ith element represents the value of x modulo the ith element of mods
       :mods   the modulo for each equation
       :modinv a function able to calculate the modular inverse of an element in R/<f>"""
    if len(eqs) != len(mods):
        raise Exception('Not a proper equation representation')
    L = len(eqs)
    for i,j in zip(range(L),range(L)):
        if i != j and gcd(mods[i].mods[j]) != mods[i].ring.one():
            raise ValueError('The moduli are not pairwise coprime')

    N = reduce((lambda x,y:x*y), mods)
    r = reduce((lambda x,y:x+y), [eqs[i]*(N/mods[i])*modinv(N/mods[i],mods[i]) for i in range(L)])
    return r
