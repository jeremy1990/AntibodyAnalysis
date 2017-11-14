from ResultLoader import loadpNovoPSMs
from ResultSorter import sortpNovoPSMs
import sys

def main():
    if len(sys.argv) != 2:
        sys.exit("Usage: python HerceptinAnalyse.py [spectracount|scorecount]")

    exprimentConditions = [
        { "enzyme": "Lys-C", "file": "H1" },
        { "enzyme": "Trypsin", "file": "H2" },
        { "enzyme": "Chymotrypsin", "file": "H3" },
        { "enzyme": "Glu-C", "file": "H4" },
        { "enzyme": "Asp-N", "file": "H5" }]
    psms = []
    for exprimentCondition in exprimentConditions:
        filePath = 'D:\\Data\\ProteomeDataAnalysis\\Antibody_H\\' \
                   + 'AnalysisResults\\' + exprimentCondition["file"] + "\\" \
                   + 'pNovoSearch\\result\\pNovo.res'
        psms.extend(loadpNovoPSMs(filePath, exprimentCondition))
    peptides = sortpNovoPSMs(psms, method=sys.argv[1], scoreRule="sum",
                             conditionKeys=["enzyme", "file"])
    printPeptides(peptides)


def printPeptides(peptides, conditionKeys=["enzyme", "file"]):
    keys = ["peptideSequence", "peptideLength", "peptideRank", "spectraCount", \
            "sumScore", "highestScore", "selectedId", "selectedTitle", \
            "selectedPrecursorMass", "selectedPrecursorError"]
    keys.extend(conditionKeys)
    for i in range(0, len(keys)-1):
        print "%s," % keys[i],
    print keys[-1]
    for peptide in peptides:
        for i in range(0, len(keys)-1):
            print "%s," % peptide[keys[i]],
        print peptide[keys[-1]]

if __name__ == "__main__":
    main()
