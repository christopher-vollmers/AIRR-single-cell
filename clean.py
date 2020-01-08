#!/usr/bin/env python3
# Roger Volden

'''
Removes the files that are completely empty
Usage:
    python3 clean.py inDir
    python3 clean.py .
'''

import sys
import os

def main():
    inDir=sys.argv[1]
    fileList = os.listdir(inDir)
    for file in fileList:
        if not os.stat(inDir+'/'+file).st_size:
            os.remove(inDir+'/'+file)

main()
