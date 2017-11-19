import unittest
import mock
import io
import antiana

class ResultLoaderTest(unittest.TestCase):
    def testLoadpNovoPSMs(self):
        filePath = 'D:\\Data\\ProteomeDataAnalysis\\Antibody_H\\AnalysisResults\\' \
                   + 'H1\\pNovoSearch\\result\\pNovo.res'
        experimentCondition = {
            "enzyme": "Lys-C"
        }
        psms = antiana.loadpNovoPSMs(filePath, experimentCondition)
        print len(psms)
