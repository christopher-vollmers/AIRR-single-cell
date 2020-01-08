#!/usr/bin/env python3
# Roger Volden

'''
This program is going to read in psl files and determine which isotype
each read in that psl aligns to based on the blocks. I'm using the
positions that I extracted from the genome browser to define where
each of the isotypes are. Each of the isotypes should have both a
membrane bound as well as a secreted form.

python3 identifyIsotypes.py isotype_positions .
'''

import sys
import os

def readIso(inFile):
    '''
    Takes my file with isotype positions and returns a dictionary
    isoDict = {'isotype':range(start, stop), ...}
    '''
    isoDict = {}
    for line in open(inFile):
        line = line.rstrip()
        if not line or line.startswith('#'): continue
        line = line.split()
        isoDict[line[0]] = [int(x) for x in line[1].split(':')[1].split('-')]
        isoDict[line[0]] = range(isoDict[line[0]][0], isoDict[line[0]][1])
    return isoDict

def readPSL(inFile):
    readDict = {}
    for line in open(inFile):
        line = line.rstrip()
        if not line: continue
        line = line.split('\t')
        header = line[9]
        bStarts = [int(x) for x in line[20].split(',')[:-1]]
        bSizes = [int(x) for x in line[18].split(',')[:-1]]

        bPairs = []
        for i in range(len(bStarts)):
            bPairs.append(set([bStarts[i], bStarts[i] + bSizes[i]]))
        readDict[header] = bPairs
    return readDict

def main():
    sub_path=sys.argv[2]
    isotypes = readIso(sys.argv[1])
    pslList = [x for x in os.listdir(sub_path) if x.endswith('-IGH.psl')]
    from progress.bar import Bar
    b = Bar('Identifying isotypes', max=len(pslList))
    for psl in pslList:
        reads = readPSL(sub_path+'/'+psl)
        outBase = psl.split('-')[0] +'-' # cell_0_barcode_
        matches = set()
        for header, pairs in reads.items():
            for pair in pairs:
                for isotype, posRange in isotypes.items():
                    if pair.intersection(posRange):
                        # print(header, isotype)
                        matches.add((header, isotype))
        for match in matches:
            out = open(sub_path+'/isotypes/' + outBase + match[1], 'a+')
            out.write(match[0] + '\t' + match[1] + '\n')
            out.close()
        b.next()
    b.finish()

main()
