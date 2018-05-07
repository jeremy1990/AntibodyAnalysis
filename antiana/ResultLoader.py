'''
This method goes through a pNovo result file and
loads all PSMs (Peptide Spectrum Matches) into an array.
Each item in the array is a PSM hash.

A typical result block which contains one PSM in pNovo is like this
S16	H1.87.87.2.dta	759.400500
P1	TKPEWV	5.518731	1.000000	-1	-4.121707	5.518731
P2	TKEPWV	5.518731	1.000000	-1	-4.121707	5.518731
P3	KTPEWV	5.518731	1.000000	-1	-4.121707	5.518731
P4	KTEPWV	5.518731	1.000000	-1	-4.121707	5.518731
P5	TKPERE	5.244681	1.000000	-3	-1.307918	5.244681
P6	TKPEER	5.244681	1.000000	-3	-1.307918	5.244681
P7	TKEPRE	5.244681	1.000000	-3	-1.307918	5.244681
P8	TKEPER	5.244681	1.000000	-3	-1.307918	5.244681
P9	KTPERE	5.244681	1.000000	-3	-1.307918	5.244681
P10	KTPEER	5.244681	1.000000	-3	-1.307918	5.244681

time:	0.003000

"S{$specNum}" line is mandantory and gives information about the
    spectrum name and the precursor mass. $num is counted from 1 as the
    result block grows. $num is not consectutively increased (mjm: bug?).
"P{$pepRank}" line is optional which gives information about the matching
    result of a peptide sequence and its related information, including
    the peptide rank ($pepRank), the peptide sequence, the PSM score,
    precursor error. This line is optional which means it could be no
    "P{$pepRank}" line or multiple lines.
"time" line is mandantory and indicates the beginning of a result block.
empty line is always existing before "time" line

After being read into memory, a PSM (Peptide Specturm Match) has below
structure:
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

Params:
    filePath the pNovo running result
    exprimentCondition the expriment condition that could be important to
                       separate results from multiple groups, such as enzyme.
                       The value should be a hash and all key-value pairs will
                       be attached to each PSM. The value is {} by default.
Return:
    An array of PSMs
'''
def loadpNovoPSMs(filePath, exprimentCondition={}):
    psms = []

    pepRank = 0
    spec = {}

    fh = open(filePath)
    try:
        lines = fh.readlines()
    finally:
        fh.close()
    for line in lines:
        line = line.rstrip('\n')
        if line[0:5] == 'time:' or line is '':
            continue
        elif line[0] == 'S':
            pepRank = 0

            parts = line.split('\t')
            if len(parts) != 3:
                raise 'Error length: %s.' % line
            spec = {
                "id": parts[0],
                "title": parts[1],
                "precursorMass": float(parts[2])
            }
        elif line[0] == 'P':
            pepRank = pepRank + 1

            parts = line.split('\t')
            if len(parts) < 7:
                raise 'Error length: %s.' % line
            if 'P' + str(pepRank) != parts[0]:
                raise Exception('Error pepRank: expected P%d, actual is %s' \
                       % (pepRank, parts[0]))
            psm = {
                "id": spec["id"] + "-" + parts[0],
                "spectrumTitle": spec["title"],
                "peptideSequence": parts[1],
                "peptideLength": len(parts[1]),
                "peptideRank": pepRank,
                "score": float(parts[2]),
                "precursorMass": spec["precursorMass"],
                "precursorError": float(parts[6])
            }
            for key in exprimentCondition:
                psm[key] = exprimentCondition[key]
            psms.append(psm)
        else:
            raise 'Error line: %s' % line
    return psms

'''
http://effbot.org/zone/default-values.htm
'''
import math
class IdentifiedProtein:
    def __init__(self, ac="", proteins=None, score=0.0, qvalue=0.0, coverage=0.0,
                       peptides=None, samesets=None, subsets=None,
                       haveDistintPep=False):
        self.ac = ac
        if proteins is None:
            self.proteins = {}
        else:
            self.proteins = proteins
        self.score = score
        self.qvalue = qvalue
        self.coverage = coverage
        if peptides is None:
            self.peptides = []
        else:
            self.peptides = peptides
        if samesets is None:
            self.samesets = {}
        else:
            self.samesets = samesets
        if subsets is None:
            self.subsets = {}
        else:
            self.subsets = subsets
        self.haveDistintPep = haveDistintPep

    def getPeptides(self, numOfMods=0):
        return [(peptide.pos, peptide.seq) for peptide in self.peptides]

    def getSequence(self):
        return self.proteins[self.ac]["seq"]

    def getMatches(self, hasHeader=True, limits=80):
        ret = ""
        
        for row in self.__range(limits):
            indexTuple = self.__index(limits, row)
            ret = ret + self.__header(limits, indexTuple) + "\n"
            ret = ret + self.getSequence()[indexTuple[0]-1:indexTuple[1]] + "\n\n"
        '''
        for peptide in self.peptides:
            ret = ret + " " * peptide.pos + peptide.seq + "\n"
        '''
        return ret

    def __range(self, limits):
        return range(0, int(math.ceil(len(self.getSequence()) * 1.0 / limits)))

    def __index(self, limits, index):
        size = len(self.__range(limits))
        remainder = len(self.getSequence()) % limits
        if index == size - 1 and remainder != 0:
            return (index * limits + 1, index * limits + remainder)
        else:
            return (index * limits + 1, (index+1) * limits)

    def __header(self, limits, indexTuple, hasSubIndex=True):
        rangeIndex = str(indexTuple[0]) + " ~ " + str(indexTuple[1])
        subIndex = ""
        for i in range(indexTuple[0]-1, indexTuple[1]):
            subIndex = subIndex + str(i % 10)
        if hasSubIndex:
            return rangeIndex + "\n" + subIndex
        else:
            return rangeIndex

'''
Sameset or Subset of an IdentifiedProtein.
'''
class SimilarProtein:
    def __init__(self, type="SameSet", ac="", coverage=0.0, peptides=None):
        self.type = type
        self.ac = ac
        self.coverage = coverage
        if peptides is None:
            self.peptides = {}
        else:
            self.peptides = peptides

class IdentifiedPeptide:
    def __init__(self, seq="", pos=-1, mods=None, finalScore=0.0, specNum=0):
        self.seq = seq
        if mods is None:
            self.mods = {}
        else:
            self.mods = mods
        self.pos = pos
        self.finalScore = finalScore
        self.specNum = specNum

def loadIdentifiedProteins(proteinFilename, databaseFilename):
    proteins = loadDatabase(databaseFilename)
    import operator

    identifiedProteins = {}

    fh = open(proteinFilename)
    try:
        lines = fh.readlines()
    finally:
        fh.close()

    numOfIdentified = 1
    identifiedProtein = None
    for line in lines:
        line = line.rstrip('\n')
        items = line.split('\t')
        if len(items) == 10 and items[0] == 'ID':
            continue
        elif len(items) == 19 and items[2] == 'ID':
            continue
        elif items[0] == str(numOfIdentified):
            numOfIdentified = numOfIdentified + 1
            identifiedProtein = IdentifiedProtein(ac = items[1], \
                proteins = proteins, score = float(items[2]), \
                qvalue = float(items[3]), coverage = float(items[4]))
            numOfPeps = int(items[5])
            numOfSameset = int(items[6])
            numOfSubset = int(items[7])
            identifiedProteins[items[1]] = identifiedProtein
        elif len(items) > 1 and items[1] == "SameSet":
            similarProtein = SimilarProtein(ac = items[2],
                coverage = float(items[3]))
            identifiedProtein.samesets[items[2]] = similarProtein
            if len(identifiedProtein.samesets) > numOfSameset:
                raise Exception('Error number of SameSet.')
        elif len(items) > 1 and items[1] == "SubSet":
            similarProtein = SimilarProtein(type = "SubSet", ac = items[2],
                coverage = float(items[3]))
            identifiedProtein.subsets[items[2]] = similarProtein
            if len(identifiedProtein.subsets) > numOfSubset:
                raise Exception('Error number of SubSet: '
                    + identifiedProtein.ac + ", "
                    + str(len(identifiedProtein.subsets)))
        elif len(items) > 2 and items[2] == str(len(identifiedProtein.peptides)+1):
            seqIL = proteins[identifiedProtein.ac]["seq"].replace('I', 'L')
            pepIL = items[3].replace('I', 'L')
            pos = seqIL.find(pepIL)
            identifiedPeptide = IdentifiedPeptide(seq = items[3], \
                specNum = int(items[26]), pos = pos)
            identifiedProtein.peptides.append(identifiedPeptide)
            if len(identifiedProtein.peptides) == numOfPeps:
                identifiedProtein.peptides = sorted(identifiedProtein.peptides, \
                    key=operator.attrgetter('pos'))
    return identifiedProteins

def loadDatabase(databaseFilename):
    proteins = {}
    ac = seq = de = None

    fh = open(databaseFilename)
    try:
        lines = fh.readlines()
    finally:
        fh.close()
    for line in lines:
        line = line.strip()
        if line[0] == '>':
            if ac is not None:
                proteins[ac] = { "seq": seq, "description": de }
            pos = line.find("\t")
            ac = line[1:pos]
            de = line[pos+1:-1]
            seq = ""
        else:
            seq = seq + line
    proteins[ac] = { "seq": seq, "description": de }

    return proteins
