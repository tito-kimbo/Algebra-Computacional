from unittest import TestSuite, TextTestRunner, defaultTestLoader
from test.z import TestZ, TestZModRing, TestZModField
import os


def testAll():
    
    test_cases = [TestZ, TestZModRing, TestZModField]

    suite = TestSuite()
    for case in test_cases:
        suite.addTests(defaultTestLoader.loadTestsFromTestCase(case))

    runner = TextTestRunner()
    runner.run(suite)

if __name__ == "__main__":
    testAll()

