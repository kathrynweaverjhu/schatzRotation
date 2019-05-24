#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.stats import norm
from statsmodels.stats.multitest import fdrcorrection
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

numSimulationIterations = 100
SVtypes = ['DEL', 'INS', 'DUP', 'INV', 'TRA']
Samples = ['N', 'Both', 'T']

observed_ensembl_files = ['/Users/cmdb/schatzRotation/SV_characteristics/annotatedBed/51Both_withTRA_regulatoryAnnotation.bed',
                          '/Users/cmdb/schatzRotation/SV_characteristics/annotatedBed/51N_only_withTRA_regulatoryAnnotation.bed',
                          '/Users/cmdb/schatzRotation/SV_characteristics/annotatedBed/51T_only_withTRA_regulatoryAnnotation.bed']
observed_vHMEC_files = ['/Users/cmdb/schatzRotation/SV_characteristics/annotatedBed/both.vHMEC_chromHMM.bed',
                        '/Users/cmdb/schatzRotation/SV_characteristics/annotatedBed/N_only.vHMEC_chromHMM.bed',
                        '/Users/cmdb/schatzRotation/SV_characteristics/annotatedBed/T_only.vHMEC_chromHMM.bed']
observed_HMEC_files = ['/Users/cmdb/schatzRotation/SV_characteristics/annotatedBed/both.HMEC_chromHMM.bed',
                       '/Users/cmdb/schatzRotation/SV_characteristics/annotatedBed/N_only.HMEC_chromHMM.bed',
                       '/Users/cmdb/schatzRotation/SV_characteristics/annotatedBed/T_only.HMEC_chromHMM.bed']
observed_myoepithelial_files = ['/Users/cmdb/schatzRotation/SV_characteristics/annotatedBed/both.myoepithelial_chromHMM.bed',
                                '/Users/cmdb/schatzRotation/SV_characteristics/annotatedBed/N_only.myoepithelial_chromHMM.bed',
                                '/Users/cmdb/schatzRotation/SV_characteristics/annotatedBed/T_only.myoepithelial_chromHMM.bed']

simulation_ensembl_files = []
simulation_vHMEC_files = []
simulation_HMEC_files = []
simulation_myoepithelial_files = []
for iterationVal in range(1,numSimulationIterations+1):
    simulation_ensembl_files.append('/Users/cmdb/schatzRotation/SV_characteristics/simulateSVs/100iterations_20190507/Iteration_{}.regulatoryAnnotation.bed'.format(str(iterationVal)))
    simulation_vHMEC_files.append('/Users/cmdb/schatzRotation/SV_characteristics/simulateSVs/100iterations_20190507/Iteration_{}.vHMEC_chromHMM.bed'.format(str(iterationVal)))
    simulation_HMEC_files.append('/Users/cmdb/schatzRotation/SV_characteristics/simulateSVs/100iterations_20190507/Iteration_{}.HMEC_chromHMM.bed'.format(str(iterationVal)))
    simulation_myoepithelial_files.append('/Users/cmdb/schatzRotation/SV_characteristics/simulateSVs/100iterations_20190507/Iteration_{}.myoepithelial_chromHMM.bed'.format(str(iterationVal)))

chromHMM_key = {'E1': 'Active TSS',
                'E2': 'Promoter Upstream TSS',
                'E3': 'Promoter Downstream TSS with DNase',
                'E4': 'Promoter Downstream TSS',
                'E5': 'Transcription 5prime',
                'E6': 'Transcription',
                'E7': 'Transcription 3prime',
                'E8': 'Weak Transcription',
                'E9': 'Transcription Regulatory',
                'E10': 'Transcription 5prime enhancer',
                'E11': 'Transcription 3prime enhancer',
                'E12': 'Transcription Weak enhancer',
                'E13': 'Active Enahncer 1',
                'E14': 'Active Enhancer 2',
                'E15': 'Active Enhancer Flank',
                'E16': 'Weak Enhancer 1',
                'E17': 'Weak Enhancer 2',
                'E18': 'Enhancer Acetylation only',
                'E19': 'DNase only',
                'E20': 'ZNF genes & repeats',
                'E21': 'Heterochromatin',
                'E22': 'Poised Promoter',
                'E23': 'Bivalent Promoter',
                'E24': 'Repressed PolyComb',
                'E25': 'Quiescent/Low'}

simplified_key = {'E1': 'promoter_flanking_region',
                  'E2': 'promoter',
                  'E3': 'promoter',
                  'E4': 'promoter',
                  'E5': 'transcription',
                  'E6': 'transcription',
                  'E7': 'transcription',
                  'E8': 'transcription',
                  'E9': 'transcription_regulatory',
                  'E10': 'enhancer',
                  'E11': 'enhancer',
                  'E12': 'enhancer',
                  'E13': 'enhancer',
                  'E14': 'enhancer',
                  'E15': 'enhancer',
                  'E16': 'enhancer',
                  'E17': 'enhancer',
                  'E18': 'enhancer',
                  'E19': 'open_chromatin_region',
                  'E20': 'heterochromatin/active_marks',
                  'E21': 'heterochromatin',
                  'E22': 'promoter',
                  'E23': 'promoter',
                  'E24': 'repressed_PolyComb',
                  'E25': 'quiescent'}
uniqueSimplified = ['promoter_flanking_region', 'promoter', 'transcription', 'transcription_regulatory', 'enhancer', 'open_chromatin_region', 'heterochromatin/active_marks','heterochromatin','repressed_PolyComb','quiescent']
uniqueEnsembl = ['CTCF_binding_site','TF_binding_site','enhancer','open_chromatin_region','promoter','promoter_flanking_region']

def dictionary_initialization(Samples, SVtypes, numIterations, RON):
    dictToReturn = {}
    for sample in Samples:
        dictToReturn[sample] = {}
        if numIterations > 0:
            for iterationVal in range(1, numIterations+1):
                dictToReturn[sample][iterationVal] = {}
                for SVtype in SVtypes:
                    if RON == 'R': #regulatory dictionary
                        dictToReturn[sample][iterationVal][SVtype] = {}
                    elif RON == 'N': #number dictionary
                        dictToReturn[sample][iterationVal][SVtype] = 0

        else:
            for SVtype in SVtypes:
                if RON == 'R': #biotype dictionary
                    dictToReturn[sample][SVtype] = {}
                elif RON == 'N': #number dictionary
                    dictToReturn[sample][SVtype] = 0
    return dictToReturn

observed_vHMEC = dictionary_initialization(Samples, SVtypes, 0, 'R')
observed_HMEC = dictionary_initialization(Samples, SVtypes, 0, 'R')
observed_myoepithelial = dictionary_initialization(Samples, SVtypes, 0, 'R')
observed_ensembl = dictionary_initialization(Samples, SVtypes, 0, 'R')
simulation_vHMEC = dictionary_initialization(Samples, SVtypes, numSimulationIterations, 'R')
simulation_HMEC = dictionary_initialization(Samples, SVtypes, numSimulationIterations, 'R')
simulation_myoepithelial = dictionary_initialization(Samples, SVtypes, numSimulationIterations, 'R')
simulation_ensembl = dictionary_initialization(Samples, SVtypes, numSimulationIterations, 'R')

def addRegulatory(dictToAddTo, sampleType, SVtype, regulatoryElement, iteration = 0):
    if iteration == 0:
        if regulatoryElement not in dictToAddTo[sampleType][SVtype]:
            dictToAddTo[sampleType][SVtype][regulatoryElement] = 0
        dictToAddTo[sampleType][SVtype][regulatoryElement] += 1

    else:
        if regulatoryElement not in dictToAddTo[sampleType][iteration][SVtype]:
            dictToAddTo[sampleType][iteration][SVtype][regulatoryElement] = 0
        dictToAddTo[sampleType][iteration][SVtype][regulatoryElement] += 1

def parseObservedAnnotationFiles(fileToParse, annotationDict, annotationType):
    fileToWriteTo = open('affected_enhancers_{}.txt'.format(annotationType), 'w+')
    sampleType = fileToParse.split('/')[-1].split('_')[0].replace('51', '')
    if 'both'.casefold() in sampleType.casefold():
        sampleType = 'Both'
    previousValues = ('', '', '', '', '')
    for i, line in enumerate(open(fileToParse),1):
        fields=line.strip('\r\n').split('\t')
        chr = fields[0]
        start = fields[1]
        SVtype = fields[3]
        ID = fields[4]
        if annotationDict == observed_ensembl:
            regulatoryElement = fields[7]
        else:
            #regulatoryElement = chromHMM_key[fields[8]]
            regulatoryElement = simplified_key[fields[8]]

        values=(chr, start, SVtype, ID, regulatoryElement)
        if previousValues != values:
            addRegulatory(annotationDict, sampleType, SVtype, regulatoryElement)
            if regulatoryElement == 'enhancer':
                fileToWriteTo.write(line.strip('\r\n') + '\n')
        previousValues = values


for fileE, fileV, fileH, fileM in zip(observed_ensembl_files, observed_vHMEC_files, observed_HMEC_files, observed_myoepithelial_files):
    parseObservedAnnotationFiles(fileE, observed_ensembl, "ensembl")
    parseObservedAnnotationFiles(fileV, observed_vHMEC, "vHMEC")
    parseObservedAnnotationFiles(fileH, observed_HMEC, "HMEC")
    parseObservedAnnotationFiles(fileM, observed_myoepithelial, "myoepithelial")

def parseSimulationAnnotationFiles(fileToParse, annotationDict):
    iterationVal = int(fileToParse.split("/")[-1].split("_")[1].split(".")[0])
    previousValues = ('', '', '','', '')
    for i, line in enumerate(open(fileToParse),1):
        fields=line.strip('\r\n').split('\t')
        chr = fields[0]
        start = fields[1]
        SVtype = fields[3]
        Sample = fields[4]
        if annotationDict == simulation_ensembl:
            regulatoryElement = fields[7]
        else:
            #regulatoryElement = chromHMM_key[fields[8]]
            regulatoryElement = simplified_key[fields[8]]

        values=(chr, start, SVtype, Sample, regulatoryElement)
        if previousValues != values:
            addRegulatory(annotationDict, Sample, SVtype, regulatoryElement, iteration = iterationVal)


for fileE, fileV, fileH, fileM in zip(simulation_ensembl_files, simulation_vHMEC_files, simulation_HMEC_files, simulation_myoepithelial_files):
    parseSimulationAnnotationFiles(fileE, simulation_ensembl)
    parseSimulationAnnotationFiles(fileV, simulation_vHMEC)
    parseSimulationAnnotationFiles(fileH, simulation_HMEC)
    parseSimulationAnnotationFiles(fileM, simulation_myoepithelial)


def getX(sample, SVtype, regulatoryElement, obs_annotationDict):
    if regulatoryElement in obs_annotationDict[sample][SVtype]:
        x = obs_annotationDict[sample][SVtype][regulatoryElement]
    else:
        x = 0
    return (x)

def computeZScore(x, u, o):
    if o == 0:
        score = XequalU(x,u)
    else:
        score = float(x-u)/o
    return(score)

def XequalU(x,u):
    if x==u:
        return True
    else:
        return False

def getU_O_Z(sample, SVtype, regulatoryElement, avgDict, stdevDict, obs_annotationDict):
    if regulatoryElement in avgDict[sample][SVtype] and regulatoryElement in stdevDict[sample][SVtype]:
        u = avgDict[sample][SVtype][regulatoryElement]
        o = stdevDict[sample][SVtype][regulatoryElement]
        z = computeZScore(getX(sample, SVtype, regulatoryElement, obs_annotationDict), u, o)
    else:
        if getX(sample, SVtype, regulatoryElement, obs_annotationDict) == 0:
            z = 'dNA'
        else:
            z = 'NA'
        u = 'NA'
        o = 'NA'
    return (u,o,z)

def append_zScore(sample, SVtype, regulatoryElement, z, dictToAppendTo, z_To_p_list1, z_To_p_list2, zAlone_list, forCorrection_list):
    dictToAppendTo[sample][SVtype][regulatoryElement] = z
    if z not in ['NA', 'dNA', True, False]:
        z_To_p_list1.append((sample, SVtype, regulatoryElement, z))
        zAlone_list.append(z)
        p_value = norm.sf(abs(z))*2
        z_To_p_list2.append((sample, SVtype, regulatoryElement, p_value))
        forCorrection_list.append(p_value)

def overall_zScore(avgDict, stdevDict, avg_intermediateDict, obs_annotationDict, Samples, SVtypes, regulatoryElementList, annotationType):
    for sample in Samples:
        for SVtype in SVtypes:
            for regulatoryElement in regulatoryElementList:
                if len(avg_intermediateDict[sample][SVtype][regulatoryElement]) > 0:
                    avgDict[sample][SVtype][regulatoryElement] = np.average(avg_intermediateDict[sample][SVtype][regulatoryElement])
                else:
                    avgDict[sample][SVtype][regulatoryElement] = 0
                if len(avg_intermediateDict[sample][SVtype][regulatoryElement]) > 1:
                    stdevDict[sample][SVtype][regulatoryElement] = np.std(avg_intermediateDict[sample][SVtype][regulatoryElement])
    zScores = {}
    zAlone = []
    z_To_p = []
    forCorrection = []
    z_To_p_2 = []
    num_zScores = 0
    num_NA_simulation = 0
    num_NA_sim_obs = 0
    num_std_True = 0
    num_std_False = 0
    file = open('rawZscores_withObsSim_regulatory_{}.txt'.format(annotationType), 'w+')
    file.write('#Sample\tSVtype\tRegulatory_Element\tAnnotation_type\tz\tObserved\tAvg_Simulation\tStd_Dev_Simulation\n')
    for sample in Samples:
        zScores[sample] = {}
        for SVtype in SVtypes:
            zScores[sample][SVtype] = {}
            for regulatoryElement in regulatoryElementList:
                x = getX(sample, SVtype, regulatoryElement, obs_annotationDict)
                u, o, z = getU_O_Z(sample, SVtype, regulatoryElement, avgDict, stdevDict, obs_annotationDict)
                append_zScore(sample, SVtype, regulatoryElement, z, zScores, z_To_p, z_To_p_2, zAlone, forCorrection)
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
                    file.write(sample + '\t' + SVtype + '\t' + regulatoryElement + '\t'+ annotationType + '\t' + str(z) + '\t' + str(x) + '\t' + str(u) + '\t' + str(o) + '\n')

    file.close()

    error_rate = 0.05
    #rejected, pvals_corrected = fdrcorrection(forCorrection, alpha=error_rate)
    #file = open('FDRcorrection{}.txt'.format(annotationType), 'w+')
    #file.write('#Sample\tSVtype\tRegulatory_Element\tCorrected_pValue_(FDR)\tNumber_Observed\n')
    rejected, pvals_corrected, alphaSidak, alphaBonf = multipletests(forCorrection, alpha = error_rate, method='bonferroni')
    file = open('bonferroniCorrection{}.txt'.format(annotationType), 'w+')
    file.write('#Sample\tSVtype\tRegulatory_Element\tCorrected_pValue_(Bonferroni)\tNumber_Observed\n')
    for rejection_value, pval_correct, tuple_value in zip(rejected, pvals_corrected, z_To_p_2):
        sample, SVtype, regulatoryElement, p_value = tuple_value
        if rejection_value == True:
            if regulatoryElement in obs_annotationDict[sample][SVtype]:
                file.write(sample + '\t' + SVtype + '\t' + regulatoryElement + '\t' + str(pval_correct) + str(obs_annotationDict[sample][SVtype][regulatoryElement]) +'\n')
            else:
                file.write(sample + '\t' + SVtype + '\t' + regulatoryElement + '\t' + str(pval_correct) + '\t0\n')

    max_value = abs(max(zAlone))
    min_value = abs(min(zAlone))

    zCurveX = np.arange(min_value, max_value, 0.001)
    axinsX = np.arange(-15, 15, 0.001)
    zCurveY = norm.pdf(zCurveX, 0, 1)
    axinsY = norm.pdf(axinsX, 0, 1)
    fig, ax = plt.subplots()
    ax.plot(zCurveX, zCurveY)
    fig.suptitle('Standard Normal Curve {}'.format(annotationType))
    ax.set_xlabel('# of Standard Deviations around the mean')
    ax.scatter(zAlone, norm.pdf(zAlone), alpha=0.4, c='black')

    axins = inset_axes(ax, 2, 2, loc=1)
    axins.set_xlim(-15,15)
    axins.set_ylim(min(axinsY), max(axinsY))
    axins.plot(axinsX, axinsY)
    axins.scatter(zAlone, norm.pdf(zAlone), alpha=0.4, c='black')

    fig.savefig('standardNormal_withZscores_{}.png'.format(annotationType))
    plt.close()

'''Boxplot of simulation iterations'''
def boxplotByRegulatory(annotationDict, annotationType, obs_annotationDict, regulatoryElementList, Samples, SVtypes):
    avg_sim_regulatory_intermediate={}
    avg_sim_regulatory={}
    stdev_sim_regulatory={}
    for sample in Samples:
        avg_sim_regulatory_intermediate[sample]={}
        avg_sim_regulatory[sample]={}
        stdev_sim_regulatory[sample]={}
        for SVtype in SVtypes:
            avg_sim_regulatory_intermediate[sample][SVtype]={}
            avg_sim_regulatory[sample][SVtype]={}
            stdev_sim_regulatory[sample][SVtype]={}
            for regulatoryElement in regulatoryElementList:
                avg_sim_regulatory_intermediate[sample][SVtype][regulatoryElement]=[]
    for sample in Samples:
        for iterationVal in range(1, numSimulationIterations+1):
            for SVtype in SVtypes:
                for regulatoryElement in regulatoryElementList:
                    if regulatoryElement in annotationDict[sample][iterationVal][SVtype]:
                        valueToAppend = annotationDict[sample][iterationVal][SVtype][regulatoryElement]
                    else:
                        valueToAppend = 0
                    avg_sim_regulatory_intermediate[sample][SVtype][regulatoryElement].append(valueToAppend)
    overall_zScore(avg_sim_regulatory, stdev_sim_regulatory, avg_sim_regulatory_intermediate, obs_annotationDict, Samples, SVtypes, regulatoryElementList, annotationType)
    for sample in Samples:
        for SVtype in SVtypes:
            for regulatoryElement in regulatoryElementList:
                if regulatoryElement not in obs_annotationDict[sample][SVtype]:
                    obs_annotationDict[sample][SVtype][regulatoryElement] = 0
    pad = 5
    valuesToPlot = []
    names = []
    for i, SVtype in enumerate(['DEL']):
    #for i, SVtype in enumerate(SVtypes):
        for j, sample in enumerate(Samples):
            for k, regulatoryElement in enumerate(regulatoryElementList):
                valuesToPlot.append(avg_sim_regulatory_intermediate[sample][SVtype][regulatoryElement])
                names.append(regulatoryElement)
    fig, (ax1, ax2, ax3) = plt.subplots(len(['DEL']), len(Samples), figsize=(40,15), sharex= True)
    #fig, axes = plt.subplots(len(SVtypes), len(Samples), figsize=(40,35), sharex=True)
    #fig.suptitle("Simulation {}".format(annotationType), fontsize=80)
    #axes[0,0].set_title(Samples[0]+' '+annotationType, size=50)
    ax1.set_title('N only', size=50)
    #axes[0,1].set_title(Samples[1]+' '+annotationType, size=50)
    ax2.set_title('Both', size=50)
    #axes[0,2].set_title(Samples[2]+' '+annotationType, size=50)
    ax3.set_title('T only', size=50)
    #for ax, SVtype in zip(axes[:,0], SVtypes):
    for ax, SVtype in zip([ax1, ax2, ax3], ['DEL']):
        ax.annotate(SVtype, xy=(0,0.5), xytext=(-ax.yaxis.labelpad-pad, 0),
                    xycoords=ax.yaxis.label, textcoords='offset points',
                    size=50, ha='right', va='center')
    for k in range(len(['DEL'])*len(Samples)):
    #for k in range(len(SVtypes)*len(Samples)):
        i = math.floor(k/len(Samples))
        j = k%len(Samples)
        valuesToPlotIndex1 = k*len(regulatoryElementList)
        valuesToPlotIndex2 = valuesToPlotIndex1 + len(regulatoryElementList)
        [ax1, ax2, ax3][j].boxplot(valuesToPlot[valuesToPlotIndex1: valuesToPlotIndex2], labels=names[valuesToPlotIndex1:valuesToPlotIndex2])
        #axes[i,j].boxplot(valuesToPlot[valuesToPlotIndex1:valuesToPlotIndex2], labels=names[valuesToPlotIndex1: valuesToPlotIndex2])
        scatter_y = []
        scatter_names = []
        for regulatoryElement in regulatoryElementList:
            scatter_names.append(regulatoryElement)
            #scatter_y.append(obs_annotationDict[Samples[j]][SVtypes[i]][regulatoryElement])
            scatter_y.append(obs_annotationDict[Samples[j]][['DEL'][i]][regulatoryElement])
        [ax1, ax2, ax3][j].scatter(np.arange(len(regulatoryElementList))+1, scatter_y, s=100)
        #axes[i,j].scatter(np.arange(len(regulatoryElementList))+1, scatter_y, s=100)
        [ax1, ax2, ax3][j].set_ylabel('Count', size=30)
        #axes[i,j].set_ylabel('Count', size=30)
        for tick in [ax1, ax2, ax3][j].xaxis.get_major_ticks():
        #for tick in axes[i,j].xaxis.get_major_ticks():
            tick.label.set_fontsize(30)
            tick.label.set_rotation(90)
        for tick in [ax1, ax2, ax3][j].yaxis.get_major_ticks():
        #for tick in axes[i,j].yaxis.get_major_ticks():
            tick.label.set_fontsize(30)
    plt.tight_layout()
    fig.savefig('boxplotbyRegulatory_{}_withObservedOverlaid_onlyDELs.png'.format(annotationType))
    #fig.savefig('boxplotbyRegulatory_{}_withObservedOverlaid.png'.format(annotationType))
    plt.close(fig)

boxplotByRegulatory(simulation_ensembl, 'ensembl', observed_ensembl, uniqueEnsembl, Samples, SVtypes)
boxplotByRegulatory(simulation_vHMEC, 'vHMEC', observed_vHMEC, uniqueSimplified, Samples, SVtypes)
boxplotByRegulatory(simulation_HMEC, 'HMEC', observed_HMEC, uniqueSimplified, Samples, SVtypes)
boxplotByRegulatory(simulation_myoepithelial, 'myoepithelial', observed_myoepithelial, uniqueSimplified, Samples, SVtypes)
