from algorithms.divisibility import gcd
from algorithms.discrete_log import discrete_log
from examples.rings import *
from utils import *
import math


"""
Calculates euler's phi function for an integer n.

:n    integer
"""
def euler_phi(n):
    if isinstance(n,int):
        n = Z.build(n)
    assuming(isinstance(n,Z.Element), "n must be an integer")
    c = 0
    for i in range(1,n.val):
        if gcd(Z.build(i),n)==1:
            c+=1
    return c

"""
Returns True if the n is prime, False otherwise.

:n    integer
"""
def is_prime_aks(n):
    if type(n) == int:
        n = Z.build(n)
    assuming(isinstance(n,Z.Element), "n must be an integer.")
    L=math.log(n.val,2)
    #is perfect power?
    for i in range(2,int(L)+1):
        if int(n.val**(1/i))**i==n.val:
            return False
    
    #find r such that o_r(n)>L^2
    r,found=2,False
    while not found:
        found = True
        for k in range(1,int(L**2)+1):
            if n**k % Z.build(r) == Z.one:
                found = False
        if not found:
            r+=1
    if n.val <= r:
        return True
    
    #is mcd(k,n)!=1
    for i in range(2,r+1):
        if gcd(Z.build(i),n) != Z.one:
            return False
    
    #test (X+a)^n = X^n+a in Z_n[X] (mod X^r-1)
    sqrt_phi_r = math.sqrt(euler_phi(r))
    R = (Z/(Z*n))["X"]
    for i in range(1,int(sqrt_phi_r*L)):
        P1 = R.build([1,i]) # (x+i)^n
        P2 = R.build([1,0])**n+R.build([i]) # x^n+i
        M = R.build([1,0])**n-R.build([1]) # x^r-1
        if P1 % M != P2:
            return False
    return True
