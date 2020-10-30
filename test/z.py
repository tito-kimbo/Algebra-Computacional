from test.structures import TestRing, TestEuclideanDomain, TestField
from examples.rings import Z

class TestZ(TestEuclideanDomain):

    ring = Z

    buildZero = [0]
    buildOne = [1]
    buildTwo = [2]

    x = Z.build(11)
    q = Z.build(2)
    d = Z.build(5)
    r = Z.one


class TestZModRing(TestRing):

    ring = Z//(6*Z)

    buildZero = [0]
    buildOne = [1]
    buildTwo = [2]


Z17 = Z//(17*Z)

class TestZModField(TestField):

    ring = Z17

    buildZero = [0]
    buildOne = [1]
    buildTwo = [2]

    x = Z17.build(11)
    q = Z17.build(3)
    d = Z17.build(15)
    r = Z17.zero
