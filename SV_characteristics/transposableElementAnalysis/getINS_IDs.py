#!/usr/bin/env python3

# supp_vec_D = {'01': 'T',
#               '10': 'N',
#               '11': 'Both'}

Tonly = open('/Users/cmdb/schatzRotation/SV_characteristics/parsingSurvivorMergedVCF/greaterThan50Lessthan50k/Tonly_ifRestrictingLen.txt')

for line in Tonly:
    fields=line.strip("\r\n").split("\t")
    chr1 = fields[0]
    start = fields[1]
    IDcol = fields[2]
    IDfields = IDcol.split("_")
    ID = IDfields[0]
    INFO = fields[7].split(";")
    label, SVtype = INFO[3].split("=")
    if SVtype == 'INS':
        end = start
    label2, SVlength = INFO[2].split("=")
    SVlength = abs(int(SVlength))
    if SVtype == 'INS' and SVlength <= 7000 and SVlength >= 5000:
        print(chr1, '\t', start, '\t', end, '\t', IDcol, '\t', ID, '\t', SVtype)
