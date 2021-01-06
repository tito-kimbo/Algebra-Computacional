from functools import reduce

from python_alcp.algorithms.divisibility import gcd, eea

def modinv(x,f):
    """Receives two elements of the same ring R and returns the inverse of x in R/<f>. <f> must be maximal.
       
       :x  element of R
       :f  element of R generating a maximal ideal
    """
    if x.ring != f.ring:
        raise ValueError('Both elements must belong to the same ring')

    g,s,t = eea(x,f)
    return s

def chinese_remainder(eqs,mods):
    """Receives as argument a list of modular equations and returns the solution. The moduli must be pairwise coprime.
       
       :eqs    the ith element represents the value of x modulo the ith element of mods
       :mods   the modulo for each equation
    """
    if len(eqs) != len(mods):
        raise Exception('Not a proper equation representation')
    L = len(eqs)
    for i,j in zip(range(L),range(L)):
        if i != j and gcd(mods[i].mods[j]) != mods[i].ring.one:
            raise ValueError('The moduli are not pairwise coprime')

    N = reduce((lambda x,y:x*y), mods)
    r = reduce((lambda x,y:x+y), [eqs[i]*(N//mods[i])*modinv(N//mods[i],mods[i]) for i in range(L)])
    return r
