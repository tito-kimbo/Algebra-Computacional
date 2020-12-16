# PyCAS

Python based Computer Algebra System.

## Introduction

Below are some examples of common operations.

### Integers

```python
>>> from examples.rings import Z
>>> Z.one
1
>>> Z.one == 1
False
>>> Z.build(128)
128
>>> i128 = Z.build(128); i32 = Z.build(32)
>>> i128+i32
160
>>> i128*i32
4096
>>> i128//i32
4
>>> 32 * Z.one == i32
True
>>> i16 = Z.build(16); i256 = Z.build(256)
>>> i16 ** 2 == i256
True
```

### Quotients

```python
>>> from examples.rings import Z
>>> Z4 = Z/(4*Z)
>>> Z4
ℤ/4ℤ
>>> q2 = Z4.build(2)
>>> q6 = Z4.build(6)
>>> q2 == q6
True
>>> q2 + q6
[8]
>>> q2 + q6 == Z4.zero
True
>>> q2 / q6
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
TypeError: unsupported operand type(s) for /: 'structures.ideals.RingQuotient.Element' and 'structures.ideals.RingQuotient.Element'
>>> # We can't divide because Z4 is not a field!
>>> Z5 = Z/(5*Z)
>>> q3 = Z5.build(3); q2_ = Z5.build(2)
>>> q2_ == q2
False
>>> print(q2_.ring)
ℤ/5ℤ
>>> Z5.one / q3 == q2_
True
```

### Polynomials

```python
>>> from examples.rings import Z
>>> ZX = Z["X"]
>>> ZX.is_euclidean()
False
>>> Z5 = Z/(Z*5)
>>> Z5X = Z5["X"]
>>> Z5X
ℤ/5ℤ[X]
>>> Z5X.is_euclidean()
True
>>> Z5X.one
([1])
>>> pol = Z5X.build([1,2,3])
>>> pol
([1] + [2]*X^1 + [3]*X^2)
>>> pol ** 2
([1] + [4]*X^1 + [0]*X^2 + [2]*X^3 + [4]*X^4)
>>> pol ** 2 // pol
([1] + [2]*X^1 + [3]*X^2)
>>> pol ** 2 // pol == pol
True
>>> from algorithms.divisibility import gcd
>>> p1 = Z5X.build([1,0,1])
>>> p2 = Z5X.build([3,2,1,3])
>>> g = gcd(p1,p2)
>>> g
([2] + [4]*X^1)
>>> p1 % g
(0)
>>> p2 % g
(0)
>>>
```

### Finite fields

```python
>>> from examples.finite_fields import FiniteField
>>> # The polynomial x^6 + x + 1 is primitive in Z2
>>> F64 = FiniteField(2, [1,1,0,0,0,0,1])
>>> alpha = F64.generator()
>>> # if you chose a primitive polynomial to build the field
>>> # then alpha is a generator
>>> alpha
[([0] + [1]*X^1)]
>>> alpha**17
[([0] + [1]*X^1 + [1]*X^2 + [1]*X^5)]
>>> # Let's test some polynomials
>>> F64Y = F64["Y"]
>>> pol = F64Y.build([F64.one, F64.zero, alpha**23, F64.one])
>>> pol
([([1])] + [([1] + [1]*X^3 + [1]*X^5)]*Y^2 + [([1])]*Y^3)
>>> # Polynomials are hard to read, because the coefficients are complex
>>> # We can change the representation of the coefficients in terms of the generator
>>> F64.setRepr("alpha")
>>> pol
(1 + alpha**23*Y^2 + 1*Y^3)
>>> F64.setRepr(None)
>>> pol
([([1])] + [([1] + [1]*X^3 + [1]*X^5)]*Y^2 + [([1])]*Y^3)
>>>
```

### Example: polynomial factorization in finite field

```python
>>> from examples.finite_fields import FiniteField
>>> from algorithms.factorization import *
>>> F9 = FiniteField(3, [2,2,1], "a")
>>> F9X = F9["x"]
>>> f = F9X.build([-F9.generator(), F9.one, F9.one, F9.one])
>>> f
2a¹ + x¹ + x² + x³
>>> # We want to find the irreducible factors of f
>>> factors = berlekamp_factorization(f)
>>> factors
[2 + (1 + 2a¹)x¹ + x², a¹ + x¹]
>>> # Let's test that this is indeed a factorization
>>> all([x.is_prime() for x in factors])
True
>>> factors[0] * factors[1]
2a¹ + x¹ + x² + x³
>>> 
```

## Developing

Read our CONTRIBUTING guide for a description of the class hierarchy and recommended practices
