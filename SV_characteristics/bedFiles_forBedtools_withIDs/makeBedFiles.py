#!/usr/bin/env python3

'''USAGE: ./makeBedFiles.py
    Needs: ../parsingSurvivorMergedVCF/Tonly_51_restrictingLen.txt
            ../parsingSurvivorMergedVCF/Nonly_51_restrictingLen.txt
            ../parsingSurvivorMergedVCF/uniqBoth_restrictingLen.txt
    Outputs: T_only_51_withTRA.bed
            N_only_51_withTRA.bed
            both_51_withTRA.bed'''


# supp_vec_D = {'01': 'T',
#               '10': 'N',
#               '11': 'Both'}

T_only_txt = open('/Users/cmdb/schatzRotation/SV_characteristics/parsingSurvivorMergedVCF/Tonly_51_restrictingLen.txt')
N_only_txt = open('/Users/cmdb/schatzRotation/SV_characteristics/parsingSurvivorMergedVCF/Nonly_51_restrictingLen.txt')
uniqueBoth_text = open('/Users/cmdb/schatzRotation/SV_characteristics/parsingSurvivorMergedVCF/uniqBoth_restrictingLen.txt')

T_only_BED = open('/Users/cmdb/schatzRotation/SV_characteristics/bedFiles_forBedtools_withIDs/T_only_51_withTRA.bed', 'w')
N_only_BED = open('/Users/cmdb/schatzRotation/SV_characteristics/bedFiles_forBedtools_withIDs/N_only_51_withTRA.bed', 'w')
Both_BED = open('/Users/cmdb/schatzRotation/SV_characteristics/bedFiles_forBedtools_withIDs/both_51_withTRA.bed', 'w')


def write_line_forBed(fields):
    write_line = ''

    chr1 = fields[0]
    write_line += chr1
    write_line += "\t"

    start = fields[1]
    ID = fields[2]

    INFO = fields[7].split(';')
    #label, supp_vec = INFO[1].split('=')
    label2, SVlength = INFO[2].split('=')
    SVlength = abs(int(SVlength))
    label2, SVtype = INFO[3].split('=')
    label3, chr2 = INFO[5].split('=')
    if SVtype == 'INS' or SVtype == 'TRA':
    #if SVtype == 'INS':
        end = start
    else:
        label4, end = INFO[6].split('=')

    if int(start) < int(end):
        write_line += start
        write_line += '\t'

        write_line += end
        write_line += '\t'
    elif int(end) < int(start):
        write_line += end
        write_line += '\t'

        write_line += start
        write_line += '\t'

    else:
        write_line += start
        write_line += '\t'

        write_line += end
        write_line += '\t'

    write_line += SVtype
    write_line += '\t'

    write_line += ID
    write_line += '\n'
    return (write_line, SVtype, SVlength)

for line in T_only_txt:
    fields=line.strip("\r\n").split("\t")
    lineToPrint, SVtype, SVlength = write_line_forBed(fields)
    #if SVtype != 'TRA' and SVtype != 'NA' and (SVlength >= 50 and SVlength <= 50000):
    if SVtype != 'NA' and ((SVlength >= 50 and SVlength <= 50000) or SVtype == 'TRA'):
        T_only_BED.write(lineToPrint)

for line in N_only_txt:
    fields=line.strip("\r\n").split("\t")
    lineToPrint, SVtype, SVlength = write_line_forBed(fields)
    #if SVtype != 'TRA' and SVtype != 'NA' and (SVlength >= 50 and SVlength <= 50000):
    if SVtype != 'NA' and ((SVlength >= 50 and SVlength <= 50000) or SVtype == 'TRA'):
        N_only_BED.write(lineToPrint)

for line in uniqueBoth_text:
    fields=line.strip("\r\n").split("\t")
    lineToPrint, SVtype, SVlength = write_line_forBed(fields)

    #if SVtype != 'TRA' and SVtype != 'NA' and (SVlength >= 50 and SVlength <= 50000):
    if SVtype != 'NA' and ((SVlength >= 50 and SVlength <= 50000) or SVtype == 'TRA'):
        Both_BED.write(lineToPrint)

for file in [T_only_BED, N_only_BED, Both_BED]:
    file.close()
