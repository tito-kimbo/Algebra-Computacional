import numpy as np

def _incustom(elem, l):
    for i in range(len(l)):
        if(l[i]==elem):
            return True,i
    return False,-1

""" Returns the discrete logarithm with base a of b using Shanks 
algorithm (baby-step giant-step)."""
def discrete_log(a,h):
    if a.ring != h.ring:
        raise ValueError('Both elements must belong to the same ring') 
    N = a.ring.order()-1 # Order of the multiplicative group
    m = int(np.ceil(np.sqrt(N)))
    # Defining a hash/ordering might make this quicker?
    alpha = [a**j for j in range(m)]
    # order of an element divides order of the (multiplicative) group
    a_inv_m = a**(N-m)
    g = h
    for i in range(m):
        found,j = _incustom(g,alpha)
        if found:
            return i*m+j
        g *= a_inv_m
    return None