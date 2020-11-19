from test.structures import TestRing, TestEuclideanDomain

class PolyTest(TestRing):
    
    buildX = None

    def testDeg(self):
        R = self.ring
        self.assertEqual(R.zero.deg(), -1)
        self.assertEqual(R.one.deg(), 0)
        self.assertEqual(R.build(*self.buildX).deg(), 1)


    def testDer(self):
        R = self.ring
        X = R.build(*self.buildX)
        self.assertEqual(X.der(), R.one)
        self.assertEqual((2*X).der(), 2*R.one)
        self.assertEqual(R.one.der(), R.zero)

    

class PolyTestED(TestEuclideanDomain, PolyTest):
    pass
