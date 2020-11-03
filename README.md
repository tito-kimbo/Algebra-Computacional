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
>>> from structures.polynomials import GetPolynomials
>>> ZX = GetPolynomials(Z)
>>> ZX.is_euclidean()
False
>>> Z5 = Z/(Z*5)
>>> Z5X = GetPolynomials(Z5)
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
>>> from examples.rings import Z
>>> from structures.polynomials import GetPolynomials
>>> Z7 = Z/(7*Z)
>>> Z7X = GetPolynomials(Z7)
>>> gen = Z7X.build([2,0,1])
>>> gen.is_prime()
True
>>> ideal = gen * Z7X
>>> ideal.is_maximal()
True
>>> F49 = Z7X / ideal
>>> F49.char()
7
>>> F49.order()
49
>>> x = F49.build([0,1])
>>> x
[([0] + [1]*X^1)]
>>> x**2
[([5])]
>>> x**3
[([0] + [5]*X^1)]
>>> min([i for i in range(1,49) if x**i == F49.one])
12
>>> # x tiene orden 12
```

## Developing

Read our CONTRIBUTING guide for a description of the class hierarchy and recommended practices
