import unittest
import antiana

class ResultSorterTest(unittest.TestCase):
    def setUp(self):
        self.psms = [
            {
                "id": "S16-P1",
                "spectrumTitle": "H1.87.87.2.dta",
                "peptideSequence": "TKPEWV",
                "peptideLength": 6,
                "peptideRank": 1,
                "score": 5.518731,
                "precursorMass": 759.400500,
                "precursorError": 5.518731
            },
            {
                "id": "S16-P2",
                "spectrumTitle": "H1.87.87.2.dta",
                "peptideSequence": "TKPEWVQK",
                "peptideLength": 8,
                "peptideRank": 1,
                "score": 90.3489234,
                "precursorMass": 1493.400500,
                "precursorError": 1.343242
            },
            {
                "id": "S17-P1",
                "spectrumTitle": "H1.88.88.2.dta",
                "peptideSequence": "TKPEWV",
                "peptideLength": 6,
                "peptideRank": 1,
                "score": 9.490204,
                "precursorMass": 759.400500,
                "precursorError": 3.423474
            }
        ]
        self.experimentCondition = {
            "enzyme": "Lys-C"
        }

    def testSortpNovoPSMs(self):
        psmsInPepView = antiana.convertToPeptideView(self.psms,
                                conditionKeys=self.experimentCondition.keys())
        peptides = antiana.sortpNovoPSMs(psmsInPepView,
                           conditionKeys=self.experimentCondition.keys())
        assert len(peptides) == 2, \
               "expected [2], actual [%d]" % len(peptides)
        assert peptides[0]["peptideSequence"] == "TKPEWV", \
               "expected [\"TKPEWV\"], actual [\"%s\"]" \
               % peptides[0]["peptideSequence"]
        assert peptides[0]["spectraCount"] == 2, \
               "expected [2], actual [%d]" % peptides[0]["spectraCount"]

    def testConvertToPeptideView(self):
        psmsInPepView = antiana.convertToPeptideView(self.psms,
                                conditionKeys=self.experimentCondition.keys())
        assert len(psmsInPepView) == 2, \
               "expected [2], actual [%d]" % len(psmsInPepView)
        assert psmsInPepView[0]["peptideSequence"] == "TKPEWV" or \
               psmsInPepView[0]["peptideSequence"] == "TKPEWVQK", \
               "expected [\"TKPEWV\"] or [\"TKPEWVQK\"], actual [\"%s\"]" \
               % psmsInPepView[0]["peptideSequence"]
