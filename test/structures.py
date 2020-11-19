from unittest import TestCase, skip
from examples.rings import Z
from structures.rings import Ring

"""
    Test all the structural methods and constructs for an implementation
    Eg: To test the field operations in Z/Z5, make a subclass of TestField
    with ring = Z/Z5 and fill the rest of the required class fields
"""

class MockRing(Ring):

    class Element(Ring.Element):
        def __add__(self, other):
            pass
        def __sub__(self, other):
            pass
        def inner_mul(self, other):
            pass
        def __eq__(self, other):
            pass
        def __str__(self):
            pass
        def __neg__(self):
            pass
        def is_unit(self):
            return True
        def is_prime(self):
            return False

    def char(self):
        return 0
    def order(self):
        return 0
    def __eq__(self, other):
        pass
    def __str__(self):
        pass
    def __init__(self):
        self.zero = self.Element()
        self.one = self.Element()


class TestRing(TestCase):

    """ 
        Usage: make a subclass for each ring you want to test
        you must provide the arguments to build elements 0, 1 and 1+1
    """

    ring = None

    buildZero = None    # args to build 0
    buildOne = None     # args to build 1
    buildTwo = None     # args to build 1+1


    def testBuildWrongArgs(self):
        self.assertRaises(TypeError, self.ring.build, MockRing.__str__)

    def testBuildIncorrect(self):
        R = self.ring
        self.assertNotEqual(R.zero, R.build(*self.buildOne))
    
    def testBuildCorrect(self):
        R = self.ring
        self.assertEqual(R.zero, R.build(*self.buildZero))
        self.assertEqual(R.one, R.build(*self.buildOne))

    def testBuildIsElement(self):
        self.assertIsInstance(self.ring.build(*self.buildTwo), self.ring.Element)

    def testElementIsInRing(self):
        R = self.ring
        self.assertEqual(R.Element.ring, R)
        self.assertEqual(R.one.ring, R)

    def testAddCorrect(self):
        R = self.ring
        self.assertEqual(R.one + R.zero, R.one)
        self.assertEqual(R.one + R.one, R.build(*self.buildTwo))
        
    def testSubCorrect(self):
        R = self.ring
        self.assertEqual(R.one - R.zero, R.one)
        self.assertEqual(R.build(*self.buildTwo) - R.one, R.one)

    def testOuterMulCorrect(self):
        R = self.ring
        two = R.build(*self.buildTwo)
        self.assertEqual(10*R.one, 5*two)
        self.assertEqual(23*R.one, 11*two+R.one)
        
    def testMulCorrect(self):
        R = self.ring
        self.assertEqual(R.one * R.zero, R.zero)
        self.assertEqual(R.build(*self.buildTwo) * R.one, R.build(*self.buildTwo))

    def testEqRandomObject(self):
        self.assertFalse(self.ring.zero == object)

    def testEqDifferentRing(self):
        self.assertFalse(self.ring.one == MockRing().one)

    def testEqDifferentElement(self):
        self.assertFalse(self.ring.one == self.ring.zero)

    def testEqCorrect(self):
        R = self.ring
        self.assertTrue(R.zero, R.zero)

    def testNegCorrect(self):
        R = self.ring
        self.assertEqual(R.zero - R.one, -R.one)

    def testStructuralEqIncorrect(self):
        self.assertFalse(self.ring == MockRing())

    def testStructuralEqCorrect(self):
        self.assertTrue(self.ring == self.ring)


class TestEuclideanDomain(TestRing):

    # Elements x, q, d, r such that x = q*d + r
    # q,r are the results of dividing x by d

    x = None
    q = None
    d = None
    r = None

    def testFloorDivCorrect(self):
        self.assertEqual(self.q, self.x // self.d)
        self.assertEqual(self.d, self.x // self.q)

    def testModCorrect(self):
        self.assertEqual(self.r, self.x % self.d)

    def testPhi(self):
        phi = self.ring.phi
        self.assertLess(phi(self.r), phi(self.x))


class TestField(TestEuclideanDomain):

    @skip("Phi is not used in fields")
    def testPhi(self):
        pass

    def testInverseCorrect(self):
        F = self.ring
        two = F.build(*self.buildTwo)
        self.assertEqual(two * two.inverse(), F.one)

    def testTrueDivCorrect(self):
        F = self.ring
        two = F.build(*self.buildTwo)
        self.assertEqual(two.inverse(), F.one / two)

