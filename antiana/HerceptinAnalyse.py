from ResultLoader import loadpNovoPSMs
from ResultLoader import loadIdentifiedProteins
from ResultSorter import sortpNovoPSMs
from ResultSorter import convertToPeptideView
import sys

def main():
    if len(sys.argv) != 2:
        sys.exit("Usage: python HerceptinAnalyse.py [spectracount|scorecount]")
    l1803()

def antibody():
    experimentConditions = [
        { "enzyme": "Lys-C", "file": "H1" },
        { "enzyme": "Trypsin", "file": "H2" },
        { "enzyme": "Chymotrypsin", "file": "H3" },
        { "enzyme": "Glu-C", "file": "H4" },
        { "enzyme": "Asp-N", "file": "H5" }]
    psms = []
    for experimentCondition in experimentConditions:
        filePath = 'D:\\Data\\ProteomeDataAnalysis\\Antibody_H\\' \
                   + 'AnalysisResults\\' + experimentCondition["file"] + "\\" \
                   + 'pNovoSearch\\result\\pNovo.res'
        psms.extend(loadpNovoPSMs(filePath, experimentCondition))
    peptides = sortpNovoPSMs(psms, method=sys.argv[1], scoreRule="sum",
                             conditionKeys=["enzyme", "file"])
    printPeptides(peptides)

def fructosyl():
    experimentConditions = [
        { "enzyme": "Trypsin", "file": "20171103_06" },
        { "enzyme": "Chymotrypsin", "file": "20171103_07" },
        { "enzyme": "Asp-N", "file": "20171103_08" },
        { "enzyme": "Glu-C", "file": "20171103_09" }
    ]
    psms = []
    for experimentCondition in experimentConditions:
        filePath = 'D:\\Data\\ProteomeDataAnalysis\\2017-11-07-AB\\' \
                   + 'AnalysisResultsFructosyl\\' + experimentCondition["file"] \
                   + '\\pNovo\\result\\pNovo.res'
        psms.extend(loadpNovoPSMs(filePath, experimentCondition))
    psmsInPepView = convertToPeptideView(psms, conditionKeys=["enzyme", "file"])
    peptides = sortpNovoPSMs(psmsInPepView, method=sys.argv[1], scoreRule="max",
                             conditionKeys=["enzyme", "file"])
    printPeptides(peptides)

def l1803():
    experimentConditions = [
        { "enzyme": "Trypsin", "file": "L1803_Y_T" },
        { "enzyme": "Trypsin", "file": "L1803_Y_T2" },
        { "enzyme": "Trypsin", "file": "L1803_Y_Trept" },
        { "enzyme": "Chymotrypsin", "file": "L1803_Y_C" },
        { "enzyme": "Glu-C", "file": "L1803_Y_G" }]
    psms = []
    for experimentCondition in experimentConditions:
        filePath = 'D:\\Data\\ProteomeDataAnalysis\\L1803\\' \
                   + experimentCondition["file"] + "\\" \
                   + 'pNovo\\' + experimentCondition["file"] + '.txt'
        psms.extend(loadpNovoPSMs(filePath, experimentCondition))
    psmsInPepView = convertToPeptideView(psms, conditionKeys=["enzyme", "file"])
    peptides = sortpNovoPSMs(psmsInPepView, method=sys.argv[1], scoreRule="max",
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
    #main()
    proteins = loadIdentifiedProteins("D:\\Data\\ProteomeDataAnalysis\\2017-11-07-AB\\" \
                 + "2018-05-06\\20171103_06\\result\\pFind.protein",
                 "D:\\Data\\ProteomeDataAnalysis\\Databases\\" \
                 + "uniprotsprot_IMGT_ALL_IL-with_answer.fasta_td.fasta")
    #print proteins['Last-time'].getSequence()
    #print proteins['Last-time'].getPeptides()
    print proteins['Last-time'].getMatches()
