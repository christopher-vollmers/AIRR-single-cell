#!/usr/bin/env python3
# Roger Volden

'''
This program is going to have to be structured a bit differently from
the other filtering script because it's hard to read in 122G worth of
subreads to query. So instead, what I'll do is collect headers from
all of the isotype files along with their isotype. From there, I'll
iterate through the subreads file and pull out the reads that I need.
I'll use the header to isotype/cell dictionary to assign it to the
correct cell and isotype.

python3 filter_subreads.py . ../../misc/41_subreads.fastq
'''

import sys
import os

def readIso(inDir):
    '''
    Takes a directory with isotype files and reads them all in.
    A dict with the header, cell, and isotype is returned.
    isoDict = {header: cell_#_bc_isotype, ...}
    '''
    isoDict = {}
    fileList = os.listdir(inDir)
    fileList = [x for x in fileList if x.endswith('_S') or x.endswith('_M')]
    for file in fileList:
        opened = open(file, 'r')
        for line in opened:
            line = line.rstrip().split()
            header = line[0].split('_')[0]
            isoDict[header] = file
        opened.close()
    return isoDict

def readSubs(inFile, isoDict):
    '''
    Takes a subread file and the isoDict to print out the reads
    into the appropriate cell and isotype files.
    '''
    lineNum = 0
    header, seq, quality = '', '', ''
    print('Starting to read in the subreads')
    for line in open(inFile):
        line = line.rstrip()
        if not line: continue
        if lineNum % 4 == 0 and line[0] == '@':
            if header:
                head = '-'.join(header.split('_')[0].split('-')[:5])
                if head in isoDict:
                    outSub = open(isoDict[head] + '_subs.fastq', 'a+')
                    outSub.write('@' + head + '\n')
                    outSub.write(seq + '\n+\n')
                    outSub.write(quality + '\n')
                    outSub.close()
                    print(isoDict[head])
            header = line[1:]
        if lineNum % 4 == 1:
            seq = line
        if lineNum % 4 == 3:
            quality = line
        lineNum += 1
    if header:
        head = '-'.join(header.split('_')[0].split('-')[:5])
        if head in isoDict:
            outSub = open(isoDict[head] + '_subs.fastq', 'a+')
            outSub.write('@' + head + '\n')
            outSub.write(seq + '\n+\n')
            outSub.write(quality + '\n')
            outSub.close()
            print(isoDict[head])

# def readsubs(infile, isodict):
#     '''
#     takes a subread file and the isoDict to print out the reads
#     into the appropriate cell and isotype files.
#     '''
#     linenum = 0
#     header, seq, quality = '', '', ''
#     print('starting to read in the subreads')
#     for line in open(inFile):
#         line = line.rstrip()
#         if not line: continue
#         if lineNum % 4 == 0 and line[0] == '@':
#             if header:
#                 if header.split('_')[0] in isoDict:
#                     outSub = open(isoDict[header.split('_')[0]] + '_subs.fastq', 'a+')
#                     outSub.write('@' + header + '\n')
#                     outSub.write(seq + '\n+\n')
#                     outSub.write(quality + '\n')
#                     outSub.close()
#                     print(isoDict[header.split('_')[0]])
#             header = line[1:]
#         if lineNum % 4 == 1:
#             seq = line
#         if lineNum % 4 == 3:
#             quality = line
#         lineNum += 1

def main():
    isoDict = readIso(sys.argv[1])
    print('Finished collecting headers')
    readSubs(sys.argv[2], isoDict)

main()
