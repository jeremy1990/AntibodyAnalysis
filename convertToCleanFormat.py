from __future__ import print_function

aaMap = {
  'Ala': 'A',
  'Arg': 'R',
  'Asn': 'N',
  'Asp': 'D',
  'Cys': 'C',
  'Glu': 'E',
  'Gln': 'Q',
  'Gly': 'G',
  'His': 'H',
  'Ile': 'I',
  'Leu': 'L',
  'Lys': 'K',
  'Met': 'M',
  'Phe': 'F',
  'Pro': 'P',
  'Ser': 'S',
  'Thr': 'T',
  'Trp': 'W',
  'Tyr': 'Y',
  'Val': 'V'
}

assert 20 == len(aaMap)

def isInAAMap(word):
  return word in aaMap

def isPureNum(word):
  return word.isdigit()

def isPrefixAA(word):
  return word[0:3] in aaMap

def isSuffixAA(word):
  return word[-3:len(word)] in aaMap

def turnOffFlag():
  global proteinFlag
  if proteinFlag is True:
    proteinFlag = False
    print('')

def turnOnFlag():
  global proteinCount
  global proteinFlag
  if proteinFlag is False:
    proteinFlag = True
    proteinCount = proteinCount + 1
    if proteinCount is not 1:
      print('\n>', end='')
    else:
      print('>', end='')

def removeDigits(word):
  return ''.join([c for c in word if not c.isdigit()])

def main():
  fh = open('input.txt')
  try:
    lines = fh.readlines()
  finally:
    fh.close()
  for line in lines:
    line = line.rstrip('\n')
    words = line.split()
    for word in words:
      if isInAAMap(word):
        turnOffFlag()
        print(aaMap[word], end='')
      elif isPureNum(word):
        turnOffFlag()
      elif isPrefixAA(word):
        turnOffFlag()
        print(aaMap[word[0:3]], end='')
      elif isSuffixAA(word):
        turnOffFlag()
        print(aaMap[word[-3:len(word)]], end='')
      else:
        turnOnFlag()
        print(removeDigits(word), end=' ')

  #print('\nTotally %d protein(s).' % proteinCount, end='')

proteinCount = 0
proteinFlag = False
if __name__ == '__main__':
    main()
