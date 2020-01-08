#!/usr/bin/env python3
# Roger Volden

'''
This is basically the same as the antibody filter, except
it's going to sort out the reads that align to the TCRA
and TCRB loci.

A: chr14:22,178,907-23,021,667
B: chr7:141,997,301-142,511,567

Usage:
    python3 tcr_filter.py input_dir output_dir
    python3 tcr_filter.py . ../TCR/
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
        outputA = outDir + '/' + file.split('.')[0] + '-TCRA.sam'
        outputB = outDir + '/' + file.split('.')[0] + '-TCRB.sam'
        # os.system('samtools sort -o {1} {0}'.format(sam, bam))
        # os.system('samtools index {0}'.format(bam))

        # TCRA
        os.system('samtools view {0} >{1} chr14:22178907-23021667'.format(bam, outputA))
        # TCRB
        os.system('samtools view {0} >{1} chr7:141997301-142511567'.format(bam, outputB))

main()
