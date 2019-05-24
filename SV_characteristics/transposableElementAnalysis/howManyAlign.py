#!/usr/bin/env python3

import argparse as ap
import fasta

TE_List = [' LINE', ' SINE', 'L1 element', 'transposon', 'transposable']

parser = ap.ArgumentParser(description='FASTA Files from eFetch')
parser.add_argument('--eFetchFiles', action='store', nargs = '+', type=str, required = True, help='Filename/path of FASTA files from eFetch')
args=parser.parse_args()
files = args.eFetchFiles

'''binary'''
num_INS_TE = 0
for file in files:
    TE_aligned = False
    reader = fasta.FASTAReader(open(file))
    while TE_aligned == False:
        for ident, sequence in reader:
            for TE_element in TE_List:
                if TE_element.casefold() in ident.casefold():
                    print(file, '\t', ident)
                    num_INS_TE += 1
                    TE_aligned = True

print('------------BREAK-------------')
print(num_INS_TE)
print('------------BREAK-------------')

'''how many out of how many'''
for file in files:
    print('-------Break------\n', file, '\n-------Break------')
    reader = fasta.FASTAReader(open(file))
    num_TE_aligned = 0
    num_total = 0
    for ident, sequence in reader:
        num_total += 1
        for TE_element in TE_List:
            if TE_element.casefold() in ident.casefold():
                num_TE_aligned += 1
                break

    print('------------BREAK-------------')
    print(file, '\t', num_TE_aligned, '\t', num_total)
    print('------------BREAK-------------')
