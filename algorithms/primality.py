from algorithms.divisibility import gcd
from algorithms.discrete_log import discrete_log
from examples.rings import *
from utils import *
from random import randint
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
        if gcd(Z.build(i),n).is_unit():
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
    

"""
Returns True if it is a possible prime and False otherwise.

:n  the input number, larger than 1
:k  the number of tests to run
"""                                                     
def is_prime_miller_rabin(n,k=500):
    if type(n) == int:
        n = Z.build(n)
    
    if n.val < 2:
        raise ValueError("Input integer must be larger than 2.")
    if n.val == 2 or n.val == 3:
        return True
    if n.val % 2 == 0:
        return False
        
    # Factor the integer as 2^r*d+1
    r,d = 0,n.val-1
    while d % 2 == 0:
        r,d = r+1,d//2
    # Run test k times
    R = Z/(Z*n)
    one,minus_one = R.build(1),R.build(n.val-1)
    for _ in range(k):
        # Miller Test
        a = R.build(randint(2,n.val-2))**d
        if a != one and a != minus_one:
            witness,i = True,0
            while witness and i<r-1:
                a,i = a*a,i+1
                if a == minus_one:
                    witness = False
            if witness:
                return False     
    return True
