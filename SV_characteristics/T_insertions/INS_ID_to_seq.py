#!/usr/bin/env python3

import subprocess

T_only_INS = '/Users/cmdb/schatzRotation/SV_characteristics/T_insertions/T_only_INS.txt'
refinedPacbioVCF = '/Users/cmdb/schatzRotation/SV_characteristics/vcfFiles/51T_Pacbio_sniffles_s11_50bp_refined_noMT_InsDel_noEND.vcf'
refinedONTVCF = '/Users/cmdb/schatzRotation/SV_characteristics/vcfFiles/51T_ONT_sniffles_s11_50bp_refined_noMT_InsDel_noEND.vcf'

Pacbio_IDs = []
ONT_IDs = []
others = []
for line in open(T_only_INS):
    if 'PACBIO' in line:
        Pacbio_IDs.append(line.strip("\r\n"))
    elif 'ONT' in line:
        ONT_IDs.append(line.strip("\r\n"))
    else:
        others.append(line.strip("\r\n"))


# for IDline in Pacbio_IDs:
#     IDfields = IDline.split('_')
#     if 'pbsv' in IDfields[0]:
#         ID = IDfields[0].split('.')[2]
#     else:
#         ID = IDfields[0]
#     COMMAND = "awk '{if ($3 == %s) {print $5}}' %s" % (ID,refinedPacbioVCF)
#     ALTsequence = subprocess.check_output(COMMAND, shell=True).decode("utf-8")
#     if len(ALTsequence)>0:
#         file = open('{}.fa'.format(IDline), 'w')
#         file.write('>{}\n{}'.format(IDline, ALTsequence))
for IDline in ONT_IDs:
    IDfields = IDline.split('_')
    if 'pbsv' in IDfields[0]:
        ID = IDfields[0].split('.')[2]
    else:
        ID = IDfields[0]
    COMMAND = "awk '{if ($3 == %s) {print $5}}' %s" %(ID, refinedONTVCF)
    #COMMAND = "awk '{if ($3 == %s) {print $5}}' %s" %(ID, refinedPacbioVCF)
    ALTsequence = subprocess.check_output(COMMAND, shell=True).decode("utf-8")
    if len(ALTsequence)>0:
        file = open('{}.fa'.format(IDline), 'w')
        file.write('>{}\n{}'.format(IDline, ALTsequence))
