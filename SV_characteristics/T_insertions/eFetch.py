#!/usr/bin/env python3

import argparse as ap
import numpy as np
from bio import Entrez
#import subprocess
#import os
Entrez.email="kweave23@jhu.edu"

parser = ap.ArgumentParser(description='Get IDs from files from blast and then search using NCBI eUtilities')
parser.add_argument('--BlastOutFiles', action='store', nargs='+', type=str, required = True, help='Filename/path of the blast out files')
args=parser.parse_args()
files = args.BlastOutFiles

for file in files:
    fileToWriteTo = open("{}.fasta".format(file), 'w+')
    IDs = []
    for line in open(file):
        fields=line.strip('\r\n').split('\t')
        ID = str(fields[1])
        IDs.append(ID)
    IDs = np.unique(IDs)

    for ID in IDs:
        handle = Entrez.efetch(db='nucleotide',id=ID ,rettype='fasta')
        fileToWriteTo.write(handle.read())
        subprocess.call(["sleep", "1"])

        #os.system("conda init bash && conda activate /home-3/kweave23@jhu.edu/miniconda2/envs/eUtilities && esearch -db nucleotide -query %s | efetch -format fasta >> %s.fasta && sleep 1 && source deactivate"%(ID, file))

        #COMMAND10 = ["esearch", "-db", "nucleotide", "-query", ID]
        #COMMAND11 = ["efetch", "-format", "fasta"]
        #COMMAND2 = ["sleep", "1"]
        #COMMAND10_out = subprocess.Popen(COMMAND10, stdout=subprocess.PIPE) #pipe this to the next subprocess.call
        #subprocess.Popen(COMMAND11, stdin = COMMAND10_out.stdout, stdout=open('{}.fasta'.format(file), 'w+')) #send this output to a file
        #subprocess.call(COMMAND2) #sleep 1 second in between queries.
