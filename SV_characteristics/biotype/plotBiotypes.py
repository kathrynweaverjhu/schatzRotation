#!/usr/bin/env python3

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import numpy as np
import math
from scipy.stats import norm
from statsmodels.stats.multitest import fdrcorrection
from statsmodels.stats.multitest import multipletests

'''Usage: ./plotBiotypes.py
    Needs: ../annotatedBed/51T_only_withTRA_geneAnnotation.bed
            ../annotatedBed/51N_only_withTRA_geneAnnotation.bed
            ../annotatedBed/51Both_withTRA_geneAnnotation.bed
            ../simulateSVs/100iterations_20190507/Iteration_#.geneAnnotation.bed'''


numSimulationIterations = 100
SVtypes = ['DEL', 'INS', 'DUP', 'INV', 'TRA']
Samples = ['N','Both','T']
geneSubsets = ['gene', 'transcript', 'CDS', 'exon', 'start_codon', 'stop_codon','five_prime_utr' , 'three_prime_utr', 'Selenocysteine']

T_only_file = '/Users/cmdb/schatzRotation/SV_characteristics/annotatedBed/longest_transcript/51T_only_withTRA_geneAnnotation.bed'
#T_only_file = '/Users/cmdb/schatzRotation/SV_characteristics/annotatedBed/51T_only_withTRA_geneAnnotation.bed'
N_only_file = '/Users/cmdb/schatzRotation/SV_characteristics/annotatedBed/longest_transcript/51N_only_withTRA_geneAnnotation.bed'
#N_only_file = '/Users/cmdb/schatzRotation/SV_characteristics/annotatedBed/51N_only_withTRA_geneAnnotation.bed'
Both_file = '/Users/cmdb/schatzRotation/SV_characteristics/annotatedBed/longest_transcript/51Both_withTRA_geneAnnotation.bed'
#Both_file = '/Users/cmdb/schatzRotation/SV_characteristics/annotatedBed/51Both_withTRA_geneAnnotation.bed'

simulation_annotation_files = []
for iterationVal in range(1,numSimulationIterations+1):
    simulation_annotation_files.append('/Users/cmdb/schatzRotation/SV_characteristics/simulateSVs/100iterations_20190507/longest_transcript/Iteration_{}.geneAnnotation.bed'.format(str(iterationVal)))
    #simulation_annotation_files.append('/Users/cmdb/schatzRotation/SV_characteristics/simulateSVs/100iterations_20190507/Iteration_{}.geneAnnotation.bed'.format(str(iterationVal)))

'''for initializing dictionaries
        BON: B - biotype dictionary
        N - number dictionary'''
def dictionary_initialization(Samples, SVtypes, numIterations, BONOG):
    dictToReturn = {}
    for sample in Samples:
        dictToReturn[sample] = {}
        if numIterations > 0:
            for iterationVal in range(1, numIterations+1):
                dictToReturn[sample][iterationVal] = {}
                for SVtype in SVtypes:
                    if BONOG == 'B': #biotype dictionary
                        dictToReturn[sample][iterationVal][SVtype] = {}
                    elif BONOG == 'N': #number dictionary
                        dictToReturn[sample][iterationVal][SVtype] = 0
                    elif BONOG == 'G': #geneID dictionary
                        dictToReturn[sample][iterationVal][SVtype] = []

        else:
            for SVtype in SVtypes:
                if BONOG == 'B': #biotype dictionary
                    dictToReturn[sample][SVtype] = {}
                elif BONOG == 'N': #number dictionary
                    dictToReturn[sample][SVtype] = 0
                elif BONOG == 'G': #geneID dictionary
                    dictToReturn[sample][SVtype] = []
    return dictToReturn

'''Initialization of dictionaries'''
observed_biotype = dictionary_initialization(Samples, SVtypes, 0, 'B')
observed_num = dictionary_initialization(Samples, SVtypes, 0, 'N')
observed_geneIDs = dictionary_initialization(Samples, SVtypes, 0, 'G')
uniq_observed_geneIDs = dictionary_initialization(Samples, SVtypes, 0, 'G')
simulation_biotype = dictionary_initialization(Samples, SVtypes, numSimulationIterations, 'B')
simulation_num = dictionary_initialization(Samples, SVtypes, numSimulationIterations, 'N')
simulation_geneIDs = dictionary_initialization(Samples, SVtypes, numSimulationIterations, 'G')
uniq_simulation_geneIDs = dictionary_initialization(Samples, SVtypes, numSimulationIterations, 'G')

'''to add a biotype to a dictionary'''
def addBiotype(dictToAddTo, sampleType, SVtype, biotype, geneSubset, iteration = 0):
    if iteration == 0:
        if biotype not in dictToAddTo[sampleType][SVtype]:
            dictToAddTo[sampleType][SVtype][biotype] = {}
        if geneSubset not in dictToAddTo[sampleType][SVtype][biotype]:
            dictToAddTo[sampleType][SVtype][biotype][geneSubset] = 0
        dictToAddTo[sampleType][SVtype][biotype][geneSubset] += 1
    else:
        if biotype not in dictToAddTo[sampleType][iteration][SVtype]:
            dictToAddTo[sampleType][iteration][SVtype][biotype] = {}
        if geneSubset not in dictToAddTo[sampleType][iteration][SVtype][biotype]:
            dictToAddTo[sampleType][iteration][SVtype][biotype][geneSubset] = 0
        dictToAddTo[sampleType][iteration][SVtype][biotype][geneSubset] += 1

'''to add a number to a dictionary counting the occurence of SVtypes'''
def addNum(dictToAddTo, sampleType, SVtype, iteration = 0):
    if iteration == 0:
        dictToAddTo[sampleType][SVtype] += 1
    else:
        dictToAddTo[sampleType][iteration][SVtype] += 1

'''to parse/add information on count to the initialized dictionaries from the Observed Annotation Files'''
def parseObservedAnnotationFiles(fileToParse, sampleType):
    previousValues = ('', '', '', '')
    previousGeneINFO = ('', '')
    SVs_geneIDs = {}
    SVs_geneIDs_exons = {}
    for i, line in enumerate(open(fileToParse),1):
        fields = line.strip("\r\n").split("\t")
        chr = fields[0]
        start = fields[1]
        SVtype = fields[3]
        ID = fields[4]

        values = (chr, start, SVtype, ID)
        if previousValues != values:
            observed_geneIDs[sampleType][SVtype].append(list(SVs_geneIDs.keys()))
            SVs_geneSubsets = geneSubsets.copy()
            SVs_geneIDs = {}
            SVs_geneIDs_exons[values] = {}

        geneSubset = fields[-7]

        geneID = ''
        biotype = ''

        INFO = fields[-1].split(";")
        for value in INFO:
            if "gene_id" in value:
                label, geneID = value.split()
                geneID = geneID.replace('"', '')
                if geneID not in SVs_geneIDs_exons[values]:
                    SVs_geneIDs_exons[values][geneID] = []
            if "gene_biotype" in value:
                label, biotype = value.split()
                biotype = biotype.replace('"','')
            if geneSubset == 'exon':
                if "exon_number" in value:
                    label, exon_number = value.split()
                    SVs_geneIDs_exons[values][geneID].append(exon_number)
            else:
                if geneID != '' and biotype != '':
                    break

        # if geneSubset == 'CDS' and biotype == 'lincRNA':
        #     print('so much yuck')
        geneINFO = (geneID, biotype)

        if previousValues == values and previousGeneINFO == geneINFO:
            if previousValues == values and geneID in SVs_geneIDs:
                if geneSubset in SVs_geneIDs[geneID]:
                    SVs_geneIDs[geneID].remove(geneSubset)
                    #print(values, geneINFO, geneSubset)
                    addBiotype(observed_biotype, sampleType, SVtype, biotype, geneSubset)
                else: #already seen that subset with that geneID with that SV
                    i += -1
                    continue
            elif previousValues == values and geneID not in SVs_geneIDs:
                SVs_geneIDs[geneID] = SVs_geneSubsets.copy()
                SVs_geneIDs[geneID].remove(geneSubset)
                #print(values, geneINFO, geneSubset)
                addBiotype(observed_biotype, sampleType, SVtype, biotype, geneSubset)
            #else of new biotype, but same geneID is unlikely and didn't occur in test run

        elif previousValues == values and previousGeneINFO != geneINFO:
            if geneID in SVs_geneIDs:
                if geneSubset in SVs_geneIDs[geneID]:
                    SVs_geneIDs[geneID].remove(geneSubset)
                    #print(values, geneINFO, geneSubset)
                    addBiotype(observed_biotype, sampleType, SVtype, biotype, geneSubset)
                else: #already seen that subset with that geneID with that SV
                    i += -1
                    continue
            else: #already seen that SV but geneID is new
                SVs_geneIDs[geneID] = SVs_geneSubsets.copy()
                SVs_geneIDs[geneID].remove(geneSubset)
                #print(values, geneINFO, geneSubset)
                addBiotype(observed_biotype, sampleType, SVtype, biotype, geneSubset)

        else: #new SV whether new gene INFO or not
            SVs_geneIDs[geneID] = SVs_geneSubsets.copy()
            SVs_geneIDs[geneID].remove(geneSubset)
            #print(values, geneINFO, geneSubset)
            addBiotype(observed_biotype, sampleType, SVtype, biotype, geneSubset)
            addNum(observed_num, sampleType, SVtype)

        previousValues = values
        previousGeneINFO = geneINFO

    #return(SVs_geneIDs_exons)

'''parsing the observed annotation files'''
parseObservedAnnotationFiles(N_only_file, 'N')
parseObservedAnnotationFiles(T_only_file, 'T')
parseObservedAnnotationFiles(Both_file, 'Both')

# observed_exons = {}
# for sample in Samples:
#     observed_exons[sample] = {}
#     for SVtype in SVtypes:
#         observed_exons[sample][SVtype]={'genes':[], 'exons_per_gene':[]}
#
# for value in N_exons:
#     chr, start, SVtype, ID = value
#     num_genes = len(N_exons[value].keys())
#     observed_exons['N'][SVtype]['genes'].append(num_genes)
#     for geneID in N_exons[value]:
#         exonsHit = N_exons[value][geneID]
#         uniq_exonsHit = np.unique(exonsHit)
#         num_exons_perGene = len(uniq_exonsHit)
#         observed_exons['N'][SVtype]['exons_per_gene'].append(num_exons_perGene)
#
# for value in T_exons:
#     chr, start, SVtype, ID = value
#     num_genes = len(T_exons[value].keys())
#     observed_exons['T'][SVtype]['genes'].append(num_genes)
#     for geneID in T_exons[value]:
#         exonsHit = T_exons[value][geneID]
#         uniq_exonsHit = np.unique(exonsHit)
#         num_exons_perGene = len(uniq_exonsHit)
#         observed_exons['T'][SVtype]['exons_per_gene'].append(num_exons_perGene)
#
# for value in Both_exons:
#     chr, start, SVtype, ID = value
#     num_genes = len(Both_exons[value].keys())
#     observed_exons['Both'][SVtype]['genes'].append(num_genes)
#     for geneID in Both_exons[value]:
#         exonsHit = Both_exons[value][geneID]
#         uniq_exonsHit = np.unique(exonsHit)
#         num_exons_perGene = len(uniq_exonsHit)
#         observed_exons['Both'][SVtype]['exons_per_gene'].append(num_exons_perGene)


'''get unique geneIDs from observed'''
# for sample in Samples:
#     for SVtype in observed_geneIDs[sample]:
#         for value in observed_geneIDs[sample][SVtype]:
#             if value not in uniq_observed_geneIDs[sample][SVtype]:
#                 uniq_observed_geneIDs[sample][SVtype].append(value)
#
# for sample in Samples:
#     for SVtype in uniq_observed_geneIDs[sample]:
#         for theList in uniq_observed_geneIDs[sample][SVtype]:
#             for value in theList:
#                 file = open('forPanther_Observed_{}_{}.txt'.format(sample, SVtype), 'a')
#                 file.write(value + '\n')
#
# file.close()

'''to parse/add information on count to the initialized dictionaries from the simulation annotation files'''
def parseSimulationAnnotationFiles(filesToParse):
    previousValues = ('', '', '', '')
    previousGeneINFO = ('', '')
    SVs_geneIDs = {}
    SVs_geneIDs_exons = {}
    for j, file in enumerate(filesToParse,1):
        SVs_geneIDs_exons[j] = {}
        for i, line in enumerate(open(file),1):
            fields = line.strip("\r\n").split("\t")
            chr = fields[0]
            start = fields[1]
            SVtype = fields[3]
            Sample = fields[4]
            if Sample not in SVs_geneIDs_exons[j]:
                SVs_geneIDs_exons[j][Sample] = {}

            values = (chr, start, SVtype, Sample)
            if previousValues != values:
                simulation_geneIDs[Sample][j][SVtype].append(list(SVs_geneIDs.keys()))
                SVs_geneSubsets = geneSubsets.copy()
                SVs_geneIDs = {}
                SVs_geneIDs_exons[j][Sample][values] = {}

            geneSubset = fields[-7]


            geneID = ''
            biotype = ''

            INFO = fields[-1].split(";")
            for value in INFO:
                if "gene_id" in value:
                    label, geneID = value.split()
                    geneID = geneID.replace('"', '')
                    if geneID not in SVs_geneIDs_exons[j][Sample][values]:
                        SVs_geneIDs_exons[j][Sample][values][geneID] = []
                if "gene_biotype" in value:
                    label, biotype = value.split()
                    biotype = biotype.replace('"','')
                if geneSubset == 'exon':
                    if "exon_number" in value:
                        label, exon_number = value.split()
                        SVs_geneIDs_exons[j][Sample][values][geneID].append(exon_number)
                else:
                    if geneID != '' and biotype != '':
                        break
            # if geneSubset == 'CDS' and biotype == 'lincRNA':
            #     print('yuck')

            geneINFO = (geneID, biotype)

            if previousValues == values and previousGeneINFO == geneINFO:
                if previousValues == values and geneID in SVs_geneIDs:
                    if geneSubset in SVs_geneIDs[geneID]:
                        SVs_geneIDs[geneID].remove(geneSubset)
                        #print(values, geneINFO, geneSubset)
                        addBiotype(simulation_biotype, Sample, SVtype, biotype, geneSubset, iteration=j)
                    else: #already seen that subset with that geneID with that SV
                        i += -1
                        continue
                elif previousValues == values and geneID not in SVs_geneIDs:
                    SVs_geneIDs[geneID] = SVs_geneSubsets.copy()
                    SVs_geneIDs[geneID].remove(geneSubset)
                    #print(values, geneINFO, geneSubset)
                    addBiotype(simulation_biotype, Sample, SVtype, biotype, geneSubset, iteration=j)
                #else of new biotype, but same geneID is unlikely and didn't occur in test run

            elif previousValues == values and previousGeneINFO != geneINFO:
                if geneID in SVs_geneIDs:
                    if geneSubset in SVs_geneIDs[geneID]:
                        SVs_geneIDs[geneID].remove(geneSubset)
                        #print(values, geneINFO, geneSubset)
                        addBiotype(simulation_biotype, Sample, SVtype, biotype, geneSubset, iteration=j)
                    else: #already seen that subset with that geneID with that SV
                        i += -1
                        continue
                else: #already seen that SV but geneID is new
                    SVs_geneIDs[geneID] = SVs_geneSubsets.copy()
                    SVs_geneIDs[geneID].remove(geneSubset)
                    #print(values, geneINFO, geneSubset)
                    addBiotype(simulation_biotype, Sample, SVtype, biotype, geneSubset, iteration=j)

            else: #new SV whether new gene INFO or not
                SVs_geneIDs[geneID] = SVs_geneSubsets.copy()
                SVs_geneIDs[geneID].remove(geneSubset)
                #print(values, geneINFO, geneSubset)
                addBiotype(simulation_biotype, Sample, SVtype, biotype, geneSubset, iteration=j)
                addNum(simulation_num, Sample, SVtype, iteration=j)

            previousValues = values
            previousGeneINFO = geneINFO

        #return(SVs_geneIDs_exons)
'''parsing the simulation annotation files'''
parseSimulationAnnotationFiles(simulation_annotation_files)

# simulation_exons_intermediate = {}
# for sample in Samples:
#     simulation_exons_intermediate[sample] = {}
#     for SVtype in SVtypes:
#         simulation_exons_intermediate[sample][SVtype] = {'genes':[], 'exons_per_gene':[]}
#
# for iterVal in simulation_exons:
#     for sample in simulation_exons[iterVal]:
#         for value in simulation_exons[iterVal][sample]:
#             chr, start, recorded_SVtype, recorded_Sample = value
#             num_genes = len(simulation_exons[iterVal][sample][value].keys())
#             simulation_exons_intermediate[sample][recorded_SVtype]['genes'].append(num_genes)
#             for i, geneID in enumerate(simulation_exons[iterVal][sample][value],1):
#                 exonsHit = simulation_exons[iterVal][sample][value][geneID]
#                 uniq_exonsHit = np.unique(exonsHit)
#                 num_exons_perGene = len(uniq_exonsHit)
#                 simulation_exons_intermediate[sample][recorded_SVtype]['exons_per_gene'].append(num_exons_perGene)

'''violin plot of num_genes&exons_per_gene'''
#first column num_genes
#second column num_exons_per_gene
#first row N
#second row T
#third row Both
#DEL_O, DEL_S, INS_O, INS_S, DUP_O, DUP_S, INV_O, INV_S, TRA_O, TRA_S
# figs = ['N', 'T', 'Both']
# columns = ['genes', 'exons_per_gene']
#
# for figName in figs:
#     fig, (ax1, ax2) = plt.subplots(1,2, figsize=(40,40))
#     fig.suptitle(figName, size=55)
#     ax1.violinplot([observed_exons[figName]['DEL']['genes'],
#                     simulation_exons_intermediate[figName]['DEL']['genes'],
#                     observed_exons[figName]['INS']['genes'],
#                     simulation_exons_intermediate[figName]['INS']['genes'],
#                     observed_exons[figName]['DUP']['genes'],
#                     simulation_exons_intermediate[figName]['DUP']['genes'],
#                     observed_exons[figName]['INV']['genes'],
#                     simulation_exons_intermediate[figName]['INV']['genes'],
#                     observed_exons[figName]['TRA']['genes'],
#                     simulation_exons_intermediate[figName]['TRA']['genes']],
#                     positions=[1,1.5, 2.5,3, 4,4.5, 5.5,6, 7,7.5],
#                     showmeans = True)
#     ax2.violinplot([observed_exons[figName]['DEL']['exons_per_gene'],
#                     simulation_exons_intermediate[figName]['DEL']['exons_per_gene'],
#                     observed_exons[figName]['INS']['exons_per_gene'],
#                     simulation_exons_intermediate[figName]['INS']['exons_per_gene'],
#                     observed_exons[figName]['DUP']['exons_per_gene'],
#                     simulation_exons_intermediate[figName]['DUP']['exons_per_gene'],
#                     observed_exons[figName]['INV']['exons_per_gene'],
#                     simulation_exons_intermediate[figName]['INV']['exons_per_gene'],
#                     observed_exons[figName]['TRA']['exons_per_gene'],
#                     simulation_exons_intermediate[figName]['TRA']['exons_per_gene']],
#                     positions=[1,1.5, 2.5,3, 4,4.5, 5.5,6, 7,7.5],
#                     showmeans = True)
#     ax1.set_ylabel("number", size=45)
#     ax2.set_ylabel("number", size=45)
#     ax1.set_xticks([1,1.5, 2.5,3, 4,4.5, 5.5,6, 7,7.5])
#     ax2.set_xticks([1,1.5, 2.5,3, 4,4.5, 5.5,6, 7,7.5])
#     ax1.tick_params('y', labelsize=45)
#     ax2.tick_params('y', labelsize=45)
#     ax1.set_xticklabels(['DEL_Obs', 'DEL_Sim', 'INS_Obs', 'INS_Sim', 'DUP_Obs', 'DUP_Sim', 'INV_Obs', 'INV_Sim', 'TRA_Obs', 'TRA_Sim'], size=45)
#     ax2.set_xticklabels(['DEL_Obs', 'DEL_Sim', 'INS_Obs', 'INS_Sim', 'DUP_Obs', 'DUP_Sim', 'INV_Obs', 'INV_Sim', 'TRA_Obs', 'TRA_Sim'], size=45)
#     for tick1, tick2 in zip(ax1.get_xticklabels(), ax2.get_xticklabels()):
#         tick1.set_rotation(90)
#         tick2.set_rotation(90)
#     ax1.set_title(columns[0], size=45)
#     ax2.set_title(columns[1], size=45)
#     plt.tight_layout()
#     fig.savefig('exons_violinPlots_{}.png'.format(figName))
#     plt.close(fig)

'''get unique geneIDs from observed'''
# for sample in Samples:
#     for iterationVal in simulation_geneIDs[sample]:
#         for SVtype in simulation_geneIDs[sample][iterationVal]:
#             for value in simulation_geneIDs[sample][iterationVal][SVtype]:
#                 if value not in uniq_simulation_geneIDs[sample][iterationVal][SVtype]:
#                     uniq_simulation_geneIDs[sample][iterationVal][SVtype].append(value)
#
# for sample in Samples:
#     for iterationVal in uniq_simulation_geneIDs[sample]:
#         for SVtype in uniq_simulation_geneIDs[sample][iterationVal]:
#             for theList in uniq_simulation_geneIDs[sample][iterationVal][SVtype]:
#                 for value in theList:
#                     file = open('forPanther_Simulation{}_{}_{}.txt'.format(iterationVal, sample, SVtype), 'a')
#                     file.write(value + '\n')
# file.close()

'''Boxplot of simulation iterations'''
biotypes = []
avg_sim_biotypes_intermediate = {}
avg_sim_biotypes = {}
stdev_sim_biotypes = {}
for sample in Samples:
    avg_sim_biotypes_intermediate[sample] = {}
    avg_sim_biotypes[sample] = {}
    stdev_sim_biotypes[sample]={}
    for iterationVal in range(1, numSimulationIterations+1):
        for SVtype in SVtypes:
            for biotype in simulation_biotype[sample][iterationVal][SVtype]:
                biotypes.append(biotype)
uBiotypes = np.unique(biotypes)

for sample in Samples:
    for SVtype in SVtypes:
        avg_sim_biotypes_intermediate[sample][SVtype] = {}
        avg_sim_biotypes[sample][SVtype] = {}
        stdev_sim_biotypes[sample][SVtype] = {}
        for biotype in uBiotypes:
            avg_sim_biotypes_intermediate[sample][SVtype][biotype] = {}
            avg_sim_biotypes[sample][SVtype][biotype] = {}
            stdev_sim_biotypes[sample][SVtype][biotype] = {}
            for geneSubset in geneSubsets:
                avg_sim_biotypes_intermediate[sample][SVtype][biotype][geneSubset] = []

for sample in Samples:
    for iterationVal in range(1, numSimulationIterations+1):
        for SVtype in SVtypes:
            for biotype in simulation_biotype[sample][iterationVal][SVtype]:
                for geneSubset in geneSubsets:
                    if geneSubset in simulation_biotype[sample][iterationVal][SVtype][biotype]:
                        valueToAppend = simulation_biotype[sample][iterationVal][SVtype][biotype][geneSubset]
                        avg_sim_biotypes_intermediate[sample][SVtype][biotype][geneSubset].append(valueToAppend)
                    else:
                        avg_sim_biotypes_intermediate[sample][SVtype][biotype][geneSubset].append(0)

'''z-score for various annotations in observed data vs simulation data
computing averages and stdev from simulation data'''
for sample in Samples:
    for SVtype in SVtypes:
        for biotype in uBiotypes:
            for geneSubset in geneSubsets:
                if len(avg_sim_biotypes_intermediate[sample][SVtype][biotype][geneSubset]) > 0:
                    avg_sim_biotypes[sample][SVtype][biotype][geneSubset] = np.average(avg_sim_biotypes_intermediate[sample][SVtype][biotype][geneSubset])
                else:
                    avg_sim_biotypes[sample][SVtype][biotype][geneSubset] = 0
                if len(avg_sim_biotypes_intermediate[sample][SVtype][biotype][geneSubset]) > 1:
                    stdev_sim_biotypes[sample][SVtype][biotype][geneSubset] = np.std(avg_sim_biotypes_intermediate[sample][SVtype][biotype][geneSubset])

'''to get x value for z-score computation'''
def getX(sample, SVtype, biotype, geneSubset):
    if biotype in observed_biotype[sample][SVtype]:
        if geneSubset in observed_biotype[sample][SVtype][biotype]:
            x = observed_biotype[sample][SVtype][biotype][geneSubset]
        else:
            x = 0
    else:
        x = 0
    return (x)

'''to compute z score'''
def computeZScore (x, u, o):
    if o == 0:
        score = XequalU(x,u)
    else:
        score = float(x-u)/o
    return(score)

'''to ask if the observed x equals the simulation mean'''
def XequalU(x, u):
    if x == u:
        return True
    else:
        return False

'''to get the simulation mean, simulation stdev, and z score (numerical or text explaining why not numerical)'''
def getU_O_Z(sample, SVtype, biotype, geneSubset):
    if geneSubset in avg_sim_biotypes[sample][SVtype][biotype] and geneSubset in stdev_sim_biotypes[sample][SVtype][biotype]:
        u = avg_sim_biotypes[sample][SVtype][biotype][geneSubset]
        o = stdev_sim_biotypes[sample][SVtype][biotype][geneSubset]
        z = computeZScore(getX(sample, SVtype, biotype, geneSubset), u, o)
    else:
        if getX(sample, SVtype, biotype, geneSubset) == 0:
            z = 'dNA'
        else:
            z = 'NA'
        u = 'NA'
        o = 'NA'
    return (u, o, z)

zScores = {}
zAlone = []
z_To_p = []
forFDR = []
z_To_p_2 = []
num_zScores = 0
num_NA_simulation = 0
num_NA_sim_obs = 0
num_std_True = 0
num_std_False = 0

'''to add the zscore to the zscore dictionary'''
def append_zScore(sample, SVtype, biotype, geneSubset, z):
    zScores[sample][SVtype][biotype][geneSubset] = z
    if z not in ['NA', 'dNA', True, False]:
        z_To_p.append((sample, SVtype, biotype, geneSubset, z))
        zAlone.append(z)
        p_value = norm.sf(abs(z))*2 #two-sided
        z_To_p_2.append((sample, SVtype, biotype, geneSubset, p_value))
        forFDR.append(p_value)

'''computing z scores and adding to zscore dictionary'''
file = open('rawZscores_withObsSim.txt', 'w')
file.write('#Sample\tSVtype\tBiotype\tGene_Subset\tz\tObserved\tAvg_Simulation\tStd_Dev_Simulation\n')
for sample in Samples:
    zScores[sample] = {}
    for SVtype in SVtypes:
        zScores[sample][SVtype] = {}
        for biotype in uBiotypes:
            zScores[sample][SVtype][biotype] = {}
            if biotype in observed_biotype[sample][SVtype]:
                for geneSubset in geneSubsets:
                    x = getX(sample, SVtype, biotype, geneSubset)
                    u, o, z = getU_O_Z(sample, SVtype, biotype, geneSubset)
                    append_zScore(sample, SVtype, biotype, geneSubset, z)
                    if z == 'NA':
                        num_NA_simulation += 1
                    elif z == 'dNA':
                        num_NA_sim_obs += 1
                    elif z == True:
                        num_std_True += 1
                    elif z == False:
                        num_std_False += 1
                    else:
                        num_zScores += 1
                        file.write(sample + '\t' + SVtype + '\t' + biotype + '\t' + geneSubset + '\t' + str(z) + '\t' + str(x) + '\t' + str(u) + '\t' + str(o) + '\n')
            else:
                for geneSubset in geneSubsets:
                    u, o, z = getU_O_Z(sample, SVtype, biotype, geneSubset)
                    append_zScore(sample, SVtype, biotype, geneSubset, z)
                    if z == 'NA':
                        num_NA_simulation += 1
                    elif z == 'dNA':
                        num_NA_sim_obs += 1
                    elif z == True:
                        num_std_True += 1
                    elif z == False:
                        num_std_False += 1
                    else:
                        num_zScores += 1
                        file.write(sample + '\t' + SVtype + '\t' + biotype + '\t' + geneSubset + '\t' + str(z) + '\t0\t' + str(u) + '\t' + str(o) + '\n')


file.close()

# print(len(z_To_p))
# print(len(z_To_p_2))
# print(len(forFDR))
# print(num_zScores)
# print(num_NA_simulation)
# print(num_NA_sim_obs)
# print(num_std_True)
# print(num_std_False)
'''plotting zscores on a standard normal curve'''
max_value = abs(max(zAlone))
min_value = abs(min(zAlone))

zCurveX = np.arange(min_value, max_value, 0.001)
axinsX = np.arange(-15,15,0.001)
zCurveY = norm.pdf(zCurveX, 0, 1)
axinsY = norm.pdf(axinsX, 0,1)
fig, ax = plt.subplots()
ax.plot(zCurveX, zCurveY)
fig.suptitle('Standard Normal Curve')
ax.set_xlabel('# of Standard Deviations around the mean')
ax.scatter(zAlone, norm.pdf(zAlone), alpha=0.4, c='black')

axins = inset_axes(ax, 2,2, loc=1)
axins.set_xlim(-15,15)
axins.set_ylim(min(axinsY), max(axinsY))
axins.plot(axinsX, axinsY)
axins.scatter(zAlone, norm.pdf(zAlone), alpha=0.4, c='black')
#mark_inset(ax, axins, loc1=2, loc2=4, fc='none', ec="0.5")


fig.savefig('standardNormal_withZscores.png')
plt.close()


'''Boxplot cont....'''
for sample in Samples:
    for SVtype in SVtypes:
        for biotype in uBiotypes:
            if biotype not in observed_biotype[sample][SVtype]:
                observed_biotype[sample][SVtype][biotype]={}
                for geneSubset in geneSubsets:
                    observed_biotype[sample][SVtype][biotype][geneSubset] = 0
            else:
                for geneSubset in geneSubsets:
                    if geneSubset not in observed_biotype[sample][SVtype][biotype]:
                        observed_biotype[sample][SVtype][biotype][geneSubset] = 0
#print(avg_sim_biotypes_intermediate['T']['INS']['lincRNA']['CDS'])
#print(observed_biotype['T']['INS']['lincRNA']['CDS'])
pad = 5 # in points
for biotype in uBiotypes:
    valuesToPlot = []
    names = []
    #for i, SVtype in enumerate(['DEL']):
    for i, SVtype in enumerate(SVtypes):
        for j, sample in enumerate(Samples):
            for k, geneSubset in enumerate(geneSubsets):
                valuesToPlot.append(avg_sim_biotypes_intermediate[sample][SVtype][biotype][geneSubset])
                names.append(geneSubset)
    #fig, (ax1, ax2, ax3) = plt.subplots(len(['DEL']), len(Samples), figsize=(40,15), sharex=True)
    fig, axes = plt.subplots(len(SVtypes), len(Samples), figsize=(40,35), sharex=True)
    #fig.suptitle("Simulation {}".format(biotype), fontsize = 80)
    axes[0,0].set_title(Samples[0] + " " + biotype, size=50)
    #ax1.set_title(Samples[0] + ' ' + biotype, size=50)
    axes[0,1].set_title(Samples[1] + " " + biotype, size=50)
    #ax2.set_title(Samples[1]+ ' ' + biotype, size=50)
    axes[0,2].set_title(Samples[2] + " " + biotype, size=50)
    #ax3.set_title(Samples[2]+ ' ' + biotype, size=50)
    #for ax, SVtype in zip([ax1, ax2, ax3], ['DEL']):
    for ax, SVtype in zip(axes[:,0], SVtypes):
        ax.annotate(SVtype, xy=(0, 0.5), xytext=(-ax.yaxis.labelpad - pad, 0),
                    xycoords=ax.yaxis.label, textcoords='offset points',
                    size=50, ha='right', va='center')
    for k in range(len(SVtypes)*len(Samples)):
    #for k in range(len(['DEL'])*len(Samples)):
        i = math.floor(k/len(Samples))
        j = k%len(Samples)
        valuesToPlotIndex1 = k*len(geneSubsets)
        valuesToPlotIndex2 = valuesToPlotIndex1 + len(geneSubsets)
        axes[i,j].boxplot(valuesToPlot[valuesToPlotIndex1:valuesToPlotIndex2], labels=names[valuesToPlotIndex1: valuesToPlotIndex2])
        #[ax1, ax2, ax3][j].boxplot(valuesToPlot[valuesToPlotIndex1:valuesToPlotIndex2], labels=names[valuesToPlotIndex1:valuesToPlotIndex2])
        scatter_y = []
        scatter_names = []
        for geneSubset in geneSubsets:
            scatter_names.append(geneSubset)
            #scatter_y.append(observed_biotype[Samples[j]][['DEL'][i]][biotype][geneSubset])
            scatter_y.append(observed_biotype[Samples[j]][SVtypes[i]][biotype][geneSubset])
        #[ax1, ax2, ax3][j].scatter([1,2,3,4,5,6,7,8,9], scatter_y, s=100)
        axes[i,j].scatter([1,2,3,4,5,6,7,8,9], scatter_y, s=100)
        # if biotype == 'lincRNA' and Samples[j]=='Both':
        #     print(valuesToPlot[valuesToPlotIndex1:valuesToPlotIndex2])
        #     print(valuesToPlotIndex1, valuesToPlotIndex2)
        #     print(valuesToPlot)
        #     print(j, i)
        #     print(scatter_y)
        #     print(scatter_names)
        #     print(geneSubsets)
        #     print(Samples[j], SVtypes[i], biotype)
        #     quit()
        #[ax1, ax2, ax3][j].set_ylabel('Count', size=30)
        axes[i,j].set_ylabel('Count', size=30)
        #for tick in [ax1, ax2, ax3][j].xaxis.get_major_ticks():
        for tick in axes[i,j].xaxis.get_major_ticks():
            tick.label.set_fontsize(30)
            tick.label.set_rotation(60)
        #for tick in [ax1, ax2, ax3][j].yaxis.get_major_ticks():
        for tick in axes[i,j].yaxis.get_major_ticks():
            tick.label.set_fontsize(30)
    plt.tight_layout()
    #fig.savefig('boxplotByBiotype_{}_withObservedOverlaid_onlyDELs.png'.format(biotype))
    fig.savefig('boxplotByBiotype_{}_withObservedOverlaid.png'.format(biotype))
    plt.close(fig)

quit()

'''multiple test corrections'''
error_rate = 0.05
#rejected, pvals_corrected = fdrcorrection(forFDR, alpha=error_rate)
#file = open('FDRcorrection.txt', 'w')
#file.write('#Sample\tSVtype\tBiotype\tGene_Subset\tCorrected_pValue_(FDR)\tNumber_Observed\n')
rejected, pvals_corrected, alphaSidak, alphaBonf = multipletests(forFDR, alpha=error_rate, method = 'bonferroni')
file = open('bonferroniCorrection.txt', 'w')
file.write('#Sample\tSVtype\tBiotype\tGene_Subset\tCorrected_pValue_(Bonferroni)\tNumber_Observed\n')
for rejection_value, pval_correct, tuple_value in zip(rejected, pvals_corrected, z_To_p_2):
    sample, SVtype, biotype, geneSubset, p_value = tuple_value
    if rejection_value == True:
        if biotype in observed_biotype[sample][SVtype] and geneSubset in observed_biotype[sample][SVtype][biotype]:
            file.write(sample + '\t' + SVtype + '\t' + biotype + '\t' + geneSubset + '\t' + str(pval_correct) + '\t' + str(observed_biotype[sample][SVtype][biotype][geneSubset])+'\n')
        else:
            file.write(sample + '\t' + SVtype + '\t' + biotype + '\t' + geneSubset + '\t' + str(pval_correct) + '\t0\n')
