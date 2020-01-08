#!/usr/bin/env python3
# Roger Volden

'''
This is basically the same as the antibody filter, except
it's going to sort out the reads that align to the IGLK
and IGLL loci.

K: chr2:89132108-90540014
L: chr22:22380156-23265691

Usage:
    python3 igl_filter.py input_dir output_dir
    python3 igl_filter.py . ../IGL/
'''

import sys
import os

def main():
    inDir = sys.argv[1]
    outDir = sys.argv[2]
    samList=[]
    for x in os.listdir(inDir):
        if x.endswith('.sam'):
            if 'cell' in x and 'merged' not in x:
                samList.append(x)
    for file in samList:
        sam = inDir+'/'+file
        bam = inDir+'/'+file.split('.')[0] + '.bam'
        outputL = outDir + '/' + file.split('.')[0] + '-IGLL.sam'
        outputK = outDir + '/' + file.split('.')[0] + '-IGLK.sam'
        print(file)
        os.system('samtools sort -o {1} {0}'.format(sam, bam))
        os.system('samtools index {0}'.format(bam))

        # IGLK
        os.system('samtools view {0} >{1} chr2:89132108-90540014'.format(bam, outputK))
        # IGLL
        os.system('samtools view {0} >{1} chr22:22380156-23265691'.format(bam, outputL))

main()
