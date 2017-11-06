'''
sortpNovoPSMs sorts given psms which is a list of PSMs dictionary objects
according to given method. This method is used for picking better matching
peptides out of pNovo search results.

Params:
    psms an array of PSMs. An example of a PSM is like This
         psm = {
             "id": "S16-P1",
             "spectrumTitle": "H1.87.87.2.dta",
             "peptideSequence": "TKPEWV",
             "peptideLength": 6,
             "peptideRank": 1,
             "score": 5.518731,
             "precursorMass": 759.400500,
             "precursorError": 5.518731
         }
    method the way to sort the array and produce the sorted result.
           Possible values are "spectracount", "scorecount". Default is
           "spectracount".
           "spectracount" sorts all peptides based on the number of spectra
           in descending order. If the number of spectra is the same then
           use the length of the peptides. If the length of peptides is still
           the same, use the score of the highest PSM. If the score of the
           highest PSM is still the same, use the alphabetic order of the
           sequence.
           "scorecount" sorts all peptides based on the score of the highest
           PSM in descending order. If the score is the same, use the number
           of spectra. If the number of spectra is the same, use the length of
           the peptide. If the length of peptide is still the same, use the
           alphabetic order of the sequence.
    conditionKeys the expriment condition that could be important to
                  separate results from multiple groups, such as enzyme.
                  The value should be a list and all key-value pairs will
                  be attached to each PSM. The list is [] by default.
Return:
    A sorted array based on the given method. Each item in the array is
    as follows
    {
        "peptideSequence": "TKPEWV",
        "peptideLength": 6,
        "peptideRank": 1,
        "spectraCount": 1,
        "highestScore": 5.518731,
        "selectedId": "S16-P1",
        "selectedTitle": "H1.87.87.2.dta",
        "selectedPrecursorMass": 759.400500,
        "selectedPrecursorError": 5.518731,
        "enzyme": "Lys-C"
    }
'''
def sortpNovoPSMs(psms, method="spectracount", conditionKeys=[]):
    peptides = convertToPeptideView(psms, conditionKeys)
    if method.lower() == "spectracount":
        peptides = multikeysort(peptides, ["-spectraCount", \
                          "-peptideLength", "-highestScore", \
                          "peptideSequence"])
    elif method.lower() == "scorecount":
        peptides = multikeysort(peptides, ["-highestScore", \
                          "-spectraCount", "-peptideLength", \
                          "peptideSequence"])
    else:
        errorMsg = ["Unrecognized method type: %s", \
                    " only support spectracount and scorecount."]
        raise Exception( ",".join(errorMsg) % method)
    return peptides

'''
This method is to sort a list of dictionary by setting multiple keys.
The idea comes from https://stackoverflow.com/questions/1143671/
python-sorting-list-of-dictionaries-by-multiple-keys.

multikeysort(peptides, ["-spectraCount", "-peptideLength", "-highestScore", \
                        "peptideSequence"])
sorted(peptides, key=lambda peptide: (\
                  -peptide["spectraCount"], -peptide["peptideLength"],
                  -peptide["highestScore"], peptide["peptideSequence"])
'''
def multikeysort(items, columns):
    from operator import itemgetter
    comparers = [((itemgetter(col[1:].strip()), -1) \
                  if col.startswith('-') \
                  else (itemgetter(col.strip()), 1)) \
                  for col in columns]
    def comparer(left, right):
        for fn, mult in comparers:
            result = cmp(fn(left), fn(right))
            if result:
                return mult * result
        else:
            return 0
    return sorted(items, cmp=comparer)

def convertToPeptideView(psms, conditionKeys):
    peptides = {}
    for psm in psms:
        if psm["peptideLength"] == 0:
            continue
        if psm["peptideSequence"] not in peptides:
            peptides[psm["peptideSequence"]] = {
                "peptideSequence": psm["peptideSequence"],
                "peptideLength": psm["peptideLength"],
                "peptideRank": 1,
                "spectraCount": 1,
                "highestScore": psm["score"],
                "selectedId": psm["id"],
                "selectedTitle": psm["spectrumTitle"],
                "selectedPrecursorMass": psm["precursorMass"],
                "selectedPrecursorError": psm["precursorError"]
            }
            for key in conditionKeys:
                if key in psm:
                    peptides[psm["peptideSequence"]][key] = psm[key]
        else:
            peptideResult = peptides[psm["peptideSequence"]]
            peptideResult["spectraCount"] = peptideResult["spectraCount"] + 1
            if psm["score"] > peptideResult["highestScore"] or \
               psm["score"] == peptideResult["highestScore"] and \
               psm["score"] < peptideResult["peptideRank"]:
                peptideResult["peptideRank"] = psm["peptideRank"]
                peptideResult["highestScore"] = psm["score"]
                peptideResult["selectedId"] = psm["id"],
                peptideResult["selectedTitle"] = psm["spectrumTitle"]
                peptideResult["selectedPrecursorMass"] = psm["precursorMass"]
                peptideResult["selectedPrecursorError"] = psm["precursorError"]
                for key in conditionKeys:
                    if key in psm:
                        peptideResult[key] = psm[key]
    return [peptide[1] for peptide in peptides.items()]

def testSortpNovoPSMs():
    from resultLoader import loadpNovoPSMs
    filePath = 'D:\\Data\\ProteomeDataAnalysis\\Antibody_H\\AnalysisResults\\' \
               + 'H1\\pNovoSearch\\result\\pNovo.res'
    exprimentCondition = {
        "enzyme": "Lys-C"
    }
    psms = loadpNovoPSMs(filePath, exprimentCondition)
    peptides = sortpNovoPSMs(psms)
    print peptides
