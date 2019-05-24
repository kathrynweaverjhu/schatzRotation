#!/usr/bin/env python3

import argparse as ap
import numpy as np
import fasta

'''Usage: ./SV_Simulation.py --SV51_summarization ../basicPlots/lengthsAndNum_restrictingLen.txt --numIterations 100
    Outputs: Iteration#.bed'''

parser = ap.ArgumentParser(description='simulate structural variants in GRCh38 using benchmarks from 51 for type/number/length of SVs')
parser.add_argument('--SV51_summarization', action='store', nargs=1, type=str, required=True, help='File with a summary of SV type, number, and length from 51')
parser.add_argument('--numIterations', action='store', nargs=1, type=int, required=False, default=[10], help='Number of iterations to run the simulation; default is 10')
args=parser.parse_args()

SV51_summary = args.SV51_summarization[0] #'/Users/cmdb/schatzRotation/SV_characteristics/basicPlots/lengthsAndNum_restrictingLen.txt'
numberOfIterations = args.numIterations[0]
SVtypes = ['DEL', 'INS', 'DUP', 'INV', 'TRA']
Samples = ['T', 'N', 'Both']
chromosomes = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10',
            '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X']

'''for dictionary initialization'''
def dictionary_initialization(Samples, SVtypes):
    dictToReturn = {}
    for SVtype in SVtypes:
        dictToReturn[SVtype] = {}
        for sample in Samples:
            dictToReturn[SVtype][sample] = []
    return dictToReturn

'''for random number generation'''
def randNumGenerator():
    randomNumber = np.random.uniform(0, 1)
    return randomNumber

'''for generation of .bed line to write to simulation file'''
def generateBedLine(chr, start, len, SVtype, Sample):
    write_line = ''
    write_line += str(chr)
    write_line += '\t'

    write_line += str(start)
    write_line += '\t'

    if SVtype == 'INS' or SVtype == 'TRA':
        end = start
    else:
        end = int(start) + abs(int(len))
        end = str(end)

    write_line += str(end)
    write_line += '\t'

    write_line += SVtype
    write_line += '\t'

    write_line += Sample
    write_line += '\n'

    return(write_line)

'''to find which chr the SV is on
    and the starting int value of that chr with respect to the concatenated genome '''
def findChr(SVloc, chromosomeBoundaries):
    chr = int(next(j for j,v in enumerate(chromosomeBoundaries) if v > SVloc) +1)
    startOfChromosome = int(chromosomeBoundaries[chr-2] + 1)
    return(chr, startOfChromosome)

'''to return whether the SV crosses the boundary of the chromosome'''
def crossBoundary(SVloc, len, chromosomeBoundaries):
    chrOfStart, startOfChromosome = findChr(SVloc, chromosomeBoundaries)
    chrOfEnd, startOfChromosome = findChr(SVloc + len, chromosomeBoundaries)
    if chrOfStart == chrOfEnd:
        return(False)
    elif chrOfStart < chrOfEnd:
        return(True)

'''to return whether we have N's in the simulated SV location'''
def nProblem(SVloc, len, GRCh38):
    if GRCh38[SVloc].casefold() == 'N'.casefold():
        return(True)
    numberOfN = GRCh38[SVloc:SVloc+len].casefold().count('N'.casefold())
    if numberOfN >= 25:
        return(True)
    else:
        return(False)

'''to simulation'''
def simulation(numberOfIterations, nums, genomeLen, GRCh38, chromosomeBoundaries):
    for k in range(numberOfIterations):
        fileToWriteTo = open('/Users/cmdb/schatzRotation/SV_characteristics/simulateSVs/Iteration{}.bed'.format(str(k+1)), 'w')
        for SVtype in nums:
            for Sample in nums[SVtype]:
                i=0
                while i < nums[SVtype][Sample]:
                    SVloc = int(randNumGenerator()*genomeLen);
                    if SVtype == 'TRA':
                        len = 0
                    else:
                        len = int(lengths[SVtype][Sample][i])
                    if nProblem(SVloc, len, GRCh38):
                        pass
                    elif crossBoundary(SVloc, len, chromosomeBoundaries):
                        pass
                    else:
                        chr, startOfChromosome = findChr(SVloc, chromosomeBoundaries)
                        if chr == 23:
                            chr = 'X'
                        if chr == 1:
                            start = SVloc
                        else:
                            start = SVloc - startOfChromosome
                        fileToWriteTo.write(generateBedLine(chr, start, len, SVtype, Sample))
                        i+=1

'''Initialize Dictionaries'''
nums = dictionary_initialization(Samples, SVtypes)
lengths = dictionary_initialization(Samples, SVtypes)

'''Parse the summary input file to population the dictionaries'''
for i, line in enumerate(open(SV51_summary)):
    fields = line.strip('\r\n').split('\t')
    if i < 12:
        SVtype = fields[1]
        Sample = fields[0].split("_")[0]

        lengths_string = fields[2]
        lengths_string2 = lengths_string.replace("[", '')
        lengths_string3 = lengths_string2.replace("]", '')
        lengths_string_finis = lengths_string3.split(',')
        for value in lengths_string_finis:
            lengths[SVtype][Sample].append(int(value))
    else:
        SVtype = fields[0].split("_")[-1]
        Sample = fields[0].split("_")[0]
        if Sample == 'Total':
            pass
        else:
            nums[SVtype][Sample] = int(fields[1])

'''concatenate the genome'''
GRCh38 = ''
chromosomeBoundaries = []
for i, chromosome in enumerate(chromosomes):
    file = '/Users/cmdb/schatzRotation/SV_characteristics/downloaded_genome/Homo_sapiens.GRCh38.dna.chromosome.{}.fa'.format(chromosome)
    reader = fasta.FASTAReader (open(file))
    for ident, sequence in reader:
        lenOfSeq = len(sequence)
        if i == 0:
            chromosomeBoundaries.append(lenOfSeq)
        else:
            chromosomeBoundaries.append((lenOfSeq+chromosomeBoundaries[-1]))
        GRCh38 += sequence
genomeLen = len(GRCh38)

'''Run the Simulation'''
simulation(numberOfIterations, nums, genomeLen, GRCh38, chromosomeBoundaries)
