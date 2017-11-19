from unittest import makeSuite, TestSuite
from .ResultLoaderTest import ResultLoaderTest
from .ResultSorterTest import ResultSorterTest

def testSuite():
    testModules = [ResultLoaderTest, ResultSorterTest]
    allSuites = []
    for module in testModules:
        allSuites.append(makeSuite(module))
    return TestSuite(tuple(allSuites))
