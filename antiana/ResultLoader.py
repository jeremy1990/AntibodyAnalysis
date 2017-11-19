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
            if len(parts) != 7:
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
