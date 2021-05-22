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
    for file1 in fileList:
        print(file1,os.stat(inDir+'/'+file1).st_size)
        if os.stat(inDir+'/'+file1).st_size==0:
            os.remove(inDir+'/'+file1)

main()
