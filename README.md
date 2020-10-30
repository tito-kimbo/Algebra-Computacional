# PyCAS

Python based Computer Algebra System.

## Introduction

Below are some examples of common operations.

### Integers

```python
>>> from examples.rings import Z
>>> Z.one
<abc.examples.rings.Z.Element object at 0x7f862aadeaf0>
>>> print(Z.one)
1
>>> print(Z.build(128))
128
>>> i128 = Z.build(128); i32 = Z.build(32)
>>> print(i128+i32)
160
>>> print(i128*i32)
4096
>>> print(i128//i32)
4
```

### Quotients

```python
>>> Z4 = Z/(4*Z)
>>> print(Z4)
ℤ/4ℤ
>>> q2 = Z4.build(2)
>>> q6 = Z4.build(6)
>>> q2 == q6
True
>>> print(q2 + q6)
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

TODO

## Developing

Read our CONTRIBUTING guide for a description of the class hierarchy and recommended practices
