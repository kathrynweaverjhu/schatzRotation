#!/usr/bin/env python3
'''USAGE: ./whatAmIWorkingWith.py
Needs: ../vcfFiles/51N_spec_51T_sens_merged.survivor1.0.6.vcf
		  ../vcfFiles/51N_sens_51T_spec_merged.survivor1.0.6.vcf
Output: Tonly_51_restrictingLen.txt
        Nonly_51_restrictingLen.txt
        nonUniqueBoth_restrictingLen.txt
        sorted_nonUniqueBoth_restrictingLen.txt
        uniqBoth_restrictingLen.txt
        T_only_INS.txt
        vennDiagram.png
the purpose of this file is to count SVs and make specific files from the survivor merged vcf files. These files can be further processed/plotted with other scripts '''

import matplotlib
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

Tsens_Nspec_VCF_File = '/Users/cmdb/schatzRotation/SV_characteristics/vcfFiles/51N_spec_51T_sens_merged.survivor1.0.6.vcf'
Tspec_Nsens_VCF_File = '/Users/cmdb/schatzRotation/SV_characteristics/vcfFiles/51N_sens_51T_spec_merged.survivor1.0.6.vcf'
#original file from Sergey
#Tspec_Nsens_VCF_File = '/Users/cmdb/schatzRotation/SV_characteristics/vcfFiles/51.survivor.51N_sens.vcf'
#files that I ran through SURVIVOR
supp_vec_D = {'01': 'T',
              '10': 'N',
              '11': 'Both'}

#original file from sergey
# supp_vec_D = {'10': 'T',
#               '01': 'N',
#               '11': 'Both'}

N_only_file = open('Nonly_51_restrictingLen.txt', 'w')
non_uniq_both_file = open('nonUniqueBoth_restrictingLen.txt', 'w')
#N_only definition and file with addition to non-unique Both
N_too_short = 0
N_count_them = 0
Both_too_short = 0
Both_count_them = 0
Tonly = 0
theOthers = 0
for line in open(Tsens_Nspec_VCF_File):
    if line.startswith("#"):
        pass
    else:
        fields=line.strip("\r\n").split("\t")
        chr1 = fields[0]
        start = int(fields[1])
        ID = fields[2]

        ALTsequence = fields[4]

        INFO = fields[7].split(';')
        N51 = fields[9].split(':')
        T51 = fields[10].split(':')

        label, supp_vec = INFO[1].split('=')
        sample = supp_vec_D[supp_vec]

        label2, SVlength = INFO[2].split('=')
        SVlength = abs(int(SVlength))
        label3, SVtype = INFO[3].split('=')
        label4, chr2 = INFO[5].split('=')

        if SVtype == 'INS':
            end = start
        else:
            label5, end = INFO[6].split('=')

        if sample == 'N' and ((SVlength >= 50 and SVlength <= 50000)or SVtype == 'TRA'):
            N_count_them += 1
            N_only_file.write(line.strip('\r\n')+'\n')
        elif sample == 'N' and not ((SVlength >= 50 and SVlength <= 50000) or SVtype == 'TRA'):
            N_too_short += 1
        elif sample == 'Both' and ((SVlength >= 50 and SVlength <= 50000) or SVtype == 'TRA'):
            Both_count_them += 1
            non_uniq_both_file.write(line.strip('\r\n')+'\n')
        elif sample == 'Both' and not ((SVlength >= 50 and SVlength <= 50000) or SVtype == "TRA"):
            Both_too_short += 1
        else:
            theOthers += 1
# print("N_ONLY: ", N_count_them)
# print("WOULD BE N_ONLY: ", N_too_short)
# print("BOTH FOR THIS ONE: ", Both_count_them)
# print("WOULD BE BOTH: ", Both_too_short)
# print("ALL OTHERS: ", theOthers)

#T_only definition and file with non-unique Both addition
T_count_them = 0
T_too_short = 0
Both_count_them = 0
Both_too_short = 0
theOthers = 0
T_only_file = open('Tonly_51_restrictingLen.txt', 'w')
file = open('T_only_INS.txt', 'w')
for line in open(Tspec_Nsens_VCF_File):
    if line.startswith("#"):
        pass
    else:
        fields=line.strip("\r\n").split("\t")
        chr1 = fields[0]
        start = int(fields[1])

        ID = fields[2]

        ALTsequence = fields[4]

        INFO = fields[7].split(';')
        N51 = fields[9].split(':')
        T51 = fields[10].split(':')

        label, supp_vec = INFO[1].split('=')
        sample = supp_vec_D[supp_vec]

        label2, SVlength = INFO[2].split('=')
        SVlength = abs(int(SVlength))
        label3, SVtype = INFO[3].split('=')
        label4, chr2 = INFO[5].split('=')

        if SVtype == 'INS':
            end = start
        else:
            label5, end = INFO[6].split('=')

        if sample == 'T' and ((SVlength >= 50 and SVlength <= 50000)or SVtype == 'TRA'):
            T_count_them += 1
            if SVtype == 'INS':
                file.write(ID + '\n')
            T_only_file.write(line.strip('\r\n')+'\n')
        elif sample == 'T' and not((SVlength >= 50 and SVlength <= 50000)or SVtype == 'TRA'):
            T_too_short += 1
        elif sample == 'Both' and ((SVlength >= 50 and SVlength <= 50000)or SVtype == 'TRA'):
            Both_count_them += 1
            non_uniq_both_file.write(line.strip('\r\n')+'\n')
        elif sample == 'Both' and not ((SVlength >= 50 and SVlength <= 50000)or SVtype == 'TRA'):
            Both_too_short += 1
        else:
            theOthers += 1
# print("T_ONLY: ", T_count_them)
# print("WOULD BE T_ONLY: ", T_too_short)
# print("BOTH FOR THIS ONE: ", Both_count_them)
# print("WOULD BE BOTH: ", Both_too_short)
# print("ALL OTHERS: ", theOthers)
for theFile in [N_only_file, T_only_file, non_uniq_both_file, file]:
    theFile.close()

'''sort non_uniq_both_file'''
sorted_non_uniq = open('sorted_nonUniqueBoth_restrictingLen.txt', 'w')
non_uniq_both_file = open('nonUniqueBoth_restrictingLen.txt')
non_uniq = non_uniq_both_file.readlines()
non_uniq.sort()
for i in range(len(non_uniq)):
    sorted_non_uniq.write(non_uniq[i])

for theFile in [sorted_non_uniq, non_uniq_both_file]:
    theFile.close()

'''make a unique file from the sorted file'''
nonUniqueBoth_text = open('sorted_nonUniqueBoth_restrictingLen.txt')
uniqBoth = open('uniqBoth_restrictingLen.txt', 'w')
previousLineID = ""
for line in nonUniqueBoth_text:
    fields=line.strip("\r\n").split("\t")
    IDfields = fields[2].split("_")
    ID = IDfields[0]
    if ID == previousLineID:
        previousLineID = ID
        pass

    else:
        previousLineID = ID
        uniqBoth.write(line.strip("\r\n")+'\n')
for theFile in [nonUniqueBoth_text, uniqBoth]:
    theFile.close()

'''venn diagram'''
font = {'family': 'Arial',
        'weight': 'bold',
        'size': 12}
matplotlib.rc('font', **font)

numN = len(open('Nonly_51_restrictingLen.txt').readlines())
numT = len(open('Tonly_51_restrictingLen.txt').readlines())
numBoth = len(open('uniqBoth_restrictingLen.txt').readlines())

fig, ax = plt.subplots()
out = venn2(subsets=(numN, numT, numBoth), set_labels = ('Normal', 'Tumor'), normalize_to=200.0)
out.get_patch_by_id('10').set_color('dodgerblue')
out.get_patch_by_id('10').set_edgecolor('none')
#out.get_patch_by_id('10').set_alpha(0.4)

out.get_patch_by_id('01').set_color('olive')
out.get_patch_by_id('01').set_edgecolor('none')
#out.get_patch_by_id('01').set_alpha(0.4)

out.get_patch_by_id('11').set_color('saddlebrown')
out.get_patch_by_id('11').set_edgecolor('none')
#out.get_patch_by_id('11').set_alpha(0.4)

lblN = out.get_label_by_id("A")
lblN.set_fontsize(20)
xN, yN = lblN.get_position()
lblN.set_position((xN-3.5, yN+16))


lblT = out.get_label_by_id("B")
lblT.set_fontsize(20)
xT, yT = lblT.get_position()
lblT.set_position((xT+2.5, yT+17))


fig.savefig('vennDiagram.png')
plt.close(fig)

print(numN)
print(numT)
print(numBoth)
