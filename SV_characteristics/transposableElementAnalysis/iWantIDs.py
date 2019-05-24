#!/usr/bin/env python3

import argparse as ap
import numpy as np

parser = ap.ArgumentParser(description='Get IDs from files from blast to use to search using NCBI eUtilities')
parser.add_argument('--BlastOutFiles', action='store', nargs = '+', type=str, required = True, help='Filename/path of the blast out files')
args=parser.parse_args()
files = args.BlastOutFiles
print(files)
print(len(files))

for file in files:
    IDs = []
    for line in open(file):
        fields = line.strip('\r\n').split('\t')
        ID = str(fields[1])
        IDs.append(ID)
    IDs=np.unique(IDs)
    numberOfIDs = len(IDs)
    print("--------------",file,": ", numberOfIDs, "--------------")

    stringToPrint = ''
    for i in range(1, len(IDs)+1):
        thisID = IDs[i-1]
        stringToPrint += thisID
        stringToPrint += ' '
    print(stringToPrint)
