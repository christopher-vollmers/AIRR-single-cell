#!/usr/bin/env python3
# Roger Volden

'''
This is going to take all of the sam files in a directory and convert
them to indexed BAM files. After indexing, all of the antibody reads
will be extracted by position (chr14:105533853-106965578).

Changing again to only take out stuff that aligns to IGHA2 and IGHM
IGHA1: chr14:105703995-105708665
IGHA2: chr14:105583731-105588395
IGHM: chr14:105851705-105856218

Usage:
    python3 antibody_filter.py input_dir output_dir
    python3 antibody_filter.py . ../antibodies/
'''

import sys
import os

def main():
    inDir = sys.argv[1]
    outDir = sys.argv[2]

    samList = [x for x in os.listdir(inDir) if x.endswith('.sam')]
    for file in samList:
        sam = inDir+'/'+file
        bam = inDir+'/'+file.split('.')[0] + '.bam'
        output = outDir + '/' + file.split('.')[0] + '-IGH.sam'
        # os.system('samtools sort -o {1} {0}'.format(sam, bam))
        # os.system('samtools index {0}'.format(bam))

        # IGH
        os.system('samtools view {0} >{1} chr14:105533853-106965578'.format(bam, output))
        # IGHA1
        # output = outDir + '/' + file.split('.')[0] + '_IGHA1.sam'
        # os.system('samtools view {0} >{1} chr14:105703995-105708665'.format(bam, output))
        # # IGHA2
        # output = outDir + '/' + file.split('.')[0] + '_IGHA2.sam'
        # os.system('samtools view {0} >{1} chr14:105583731-105588395'.format(bam, output))
        # # IGHM
        # output = outDir + '/' + file.split('.')[0] + '_IGHM.sam'
        # os.system('samtools view {0} >{1} chr14:105851705-105856218'.format(bam, output))
        # # IGHA2 and IGHM
        # output = outDir + '/' + file.split('.')[0] + '_both.sam'
        # os.system('samtools view {0} >{1} chr14:105583731-105856218'.format(bam, output))

main()
