#!/usr/bin/env python3
# Roger Volden

'''
Runs sam2psl on every sam in a directory
python3 cell_sam2psl.py /media/vulpter/f3aa97db-ba06-4df6-86ba-65eb87050df1/10X_PBMCs/Promethion/4rep1/retry_mm_alignments
'''

import sys, os

inDir = sys.argv[1]
fileList = os.listdir(inDir)
fileList = [x for x in fileList if x.endswith('.sam')]

for file in fileList:
    out = inDir+'/'+file.split('.')[0] + '.psl'
    print(out)
    os.system('emtrey -m -i {0} >{1}'.format(inDir+'/'+file, out))
