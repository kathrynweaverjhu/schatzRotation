#!/usr/bin/env python3

'''Usage: ./basicPlots51vcf.py
    Needs: ../parsingSurvivorMergedVCF/uniqBoth_restrictingLen.txt
            ../parsingSurvivorMergedVCF/Tonly_51_restrictingLen.txt
            ../parsingSurvivorMergedVCF/Nonly_51_restrictingLen.txt
    Output: a bunch of png's for number (relative and raw) as well as lengths
            legnthsAndNum_restrictingLen.txt
    '''


import matplotlib
import matplotlib.pyplot as plt
import numpy as np

uniqueBothFile = '/Users/cmdb/schatzRotation/SV_characteristics/parsingSurvivorMergedVCF/uniqBoth_restrictingLen.txt'
TonlyFile = '/Users/cmdb/schatzRotation/SV_characteristics/parsingSurvivorMergedVCF/Tonly_51_restrictingLen.txt'
NonlyFile = '/Users/cmdb/schatzRotation/SV_characteristics/parsingSurvivorMergedVCF/Nonly_51_restrictingLen.txt'

def parseTheLine(fields):
    chr1 = fields[0]
    start = int(fields[1])

    INFO = fields[7].split(';')
    label, supp_vec = INFO[1].split('=')
    label2, SVlength = INFO[2].split('=')
    SVlength = abs(int(SVlength))
    label3, SVtype = INFO[3].split('=')
    label4, chr2 = INFO[5].split('=')
    if SVtype == 'INS':
        end = start
    else:
        label5, end = INFO[6].split('=')

    return(chr1, start, chr2, end, SVtype, SVlength)

T = {'chrom1':[], 'pos1':[], 'chrom2':[], 'pos2':[], 'SVtype':[], 'SVlength':[]}
for line in open(TonlyFile):
    fields=line.strip("\r\n").split("\t")
    chr1, start, chr2, end, SVtype, SVlength = parseTheLine(fields)

    if (SVlength >= 50 and SVlength <= 50000)  or SVtype =='TRA':
        T['chrom1'].append(chr1)
        T['pos1'].append(start)
        T['chrom2'].append(chr2)
        T['pos2'].append(int(end))
        T['SVtype'].append(SVtype)
        T['SVlength'].append(SVlength)

N = {'chrom1':[], 'pos1':[], 'chrom2':[], 'pos2':[], 'SVtype':[], 'SVlength':[]}
for line in open(NonlyFile):
    fields=line.strip("\r\n").split("\t")
    chr1, start, chr2, end, SVtype, SVlength = parseTheLine(fields)

    if (SVlength >= 50 and SVlength <= 50000) or SVtype =='TRA':
        N['chrom1'].append(chr1)
        N['pos1'].append(start)
        N['chrom2'].append(chr2)
        N['pos2'].append(int(end))
        N['SVtype'].append(SVtype)
        N['SVlength'].append(SVlength)

Both = {'chrom1':[], 'pos1':[], 'chrom2':[], 'pos2':[], 'SVtype':[], 'SVlength':[]}
for line in open(uniqueBothFile):
    fields=line.strip("\r\n").split("\t")
    chr1, start, chr2, end, SVtype, SVlength = parseTheLine(fields)

    if (SVlength >= 50 and SVlength <= 50000) or SVtype =='TRA':
        Both['chrom1'].append(chr1)
        Both['pos1'].append(start)
        Both['chrom2'].append(chr2)
        Both['pos2'].append(int(end))
        Both['SVtype'].append(SVtype)
        Both['SVlength'].append(SVlength)

T_type_length={'DEL':[], 'DUP':[], 'INV':[], 'TRA':[], 'INS':[], 'NA':[]}
T_lengths=[]
for type, length in zip(T['SVtype'], T['SVlength']):
        T_type_length[type].append(int(length))
        if type != 'TRA' and type != 'NA':
            T_lengths.append(length)
print("T_SVtype_and_Lengths\tDEL\t", T_type_length['DEL'])
print("T_SVtype_and_Lengths\tDUP\t", T_type_length['DUP'])
print("T_SVtype_and_Lengths\tINV\t", T_type_length['INV'])
print("T_SVtype_and_Lengths\tINS\t", T_type_length['INS'])

N_type_length={'DEL':[], 'DUP':[], 'INV':[], 'TRA':[], 'INS':[], 'NA':[]}
N_lengths=[]
for type, length in zip(N['SVtype'], N['SVlength']):
        N_type_length[type].append(int(length))
        if type != 'TRA' and type != 'NA':
            N_lengths.append(length)
print("N_SVtype_and_Lengths\tDEL\t", N_type_length['DEL'])
print("N_SVtype_and_Lengths\tDUP\t", N_type_length['DUP'])
print("N_SVtype_and_Lengths\tINV\t", N_type_length['INV'])
print("N_SVtype_and_Lengths\tINS\t", N_type_length['INS'])

Both_type_length={'DEL':[], 'DUP':[], 'INV':[], 'TRA':[], 'INS':[], 'NA':[]}
Both_lengths = []
for type, length in zip(Both['SVtype'], Both['SVlength']):
        Both_type_length[type].append(int(length))
        if type != 'TRA' and type != 'NA':
            Both_lengths.append(length)
print("Both_SVtype_and_Lengths\tDEL\t", Both_type_length['DEL'])
print("Both_SVtype_and_Lengths\tDUP\t", Both_type_length['DUP'])
print("Both_SVtype_and_Lengths\tINV\t", Both_type_length['INV'])
print("Both_SVtype_and_Lengths\tINS\t", Both_type_length['INS'])

N_num_Del = len(N_type_length['DEL'])
print("N_DEL\t", N_num_Del)
N_num_Dup = len(N_type_length['DUP'])
print("N_DUP\t", N_num_Dup)
N_num_Inv = len(N_type_length['INV'])
print("N_INV\t", N_num_Inv)
N_num_Tra = len(N_type_length['TRA'])
print("N_TRA\t", N_num_Tra)
N_num_Ins = len(N_type_length['INS'])
print("N_INS\t", N_num_Ins)
sumN = N_num_Del + N_num_Dup + N_num_Inv + N_num_Tra + N_num_Ins
N_rel_Del = N_num_Del/sumN
N_rel_Dup = N_num_Dup/sumN
N_rel_Inv = N_num_Inv/sumN
N_rel_Tra = N_num_Tra/sumN
N_rel_Ins = N_num_Ins/sumN

T_num_Del = len(T_type_length['DEL'])
print("T_DEL\t", T_num_Del)
T_num_Dup = len(T_type_length['DUP'])
print("T_DUP\t", T_num_Dup)
T_num_Inv = len(T_type_length['INV'])
print("T_INV\t", T_num_Inv)
T_num_Tra = len(T_type_length['TRA'])
print("T_TRA\t", T_num_Tra)
T_num_Ins = len(T_type_length['INS'])
print("T_INS\t", T_num_Ins)
sumT = T_num_Del + T_num_Dup + T_num_Inv + T_num_Tra + T_num_Ins
T_rel_Del = T_num_Del/sumT
T_rel_Dup = T_num_Dup/sumT
T_rel_Inv = T_num_Inv/sumT
T_rel_Tra = T_num_Tra/sumT
T_rel_Ins = T_num_Ins/sumT

Both_num_Del = len(Both_type_length['DEL'])
print("Both_DEL\t", Both_num_Del)
Both_num_Dup = len(Both_type_length['DUP'])
print("Both_DUP\t", Both_num_Dup)
Both_num_Inv = len(Both_type_length['INV'])
print("Both_INV\t", Both_num_Inv)
Both_num_Tra = len(Both_type_length['TRA'])
print("Both_TRA\t", Both_num_Tra)
Both_num_Ins = len(Both_type_length['INS'])
print("Both_INS\t", Both_num_Ins)

sumBoth = Both_num_Del + Both_num_Dup + Both_num_Inv + Both_num_Tra + Both_num_Ins
Both_rel_Del = Both_num_Del/sumBoth
Both_rel_Dup = Both_num_Dup/sumBoth
Both_rel_Inv = Both_num_Inv/sumBoth
Both_rel_Tra = Both_num_Tra/sumBoth
Both_rel_Ins = Both_num_Ins/sumBoth

print("Total_DEL\t", (N_num_Del+T_num_Del+Both_num_Del))
print("Total_DUP\t", (N_num_Dup+T_num_Dup+Both_num_Dup))
print("Total_INV\t", (N_num_Inv+T_num_Inv+Both_num_Inv))
print("Total_TRA\t", (N_num_Tra+T_num_Tra+Both_num_Tra))
print("Total_INS\t", (N_num_Ins+T_num_Ins+Both_num_Ins))

# '''vertical stacked barplot of the count of each sv type in sample type'''
# fig, ax = plt.subplots()
# Deltier = np.array((N_num_Del, T_num_Del, Both_num_Del))
# Duptier = np.array((N_num_Dup, T_num_Dup, Both_num_Dup))
# Invtier = np.array((N_num_Inv, T_num_Inv, Both_num_Inv))
# Tratier = np.array((N_num_Tra, T_num_Tra, Both_num_Tra))
# Instier = np.array((N_num_Ins, T_num_Ins, Both_num_Ins))
# plt1 = plt.bar([1,2,3], Tratier, bottom=[0,0,0], color='darkturquoise', tick_label = ["Normal only", "Tumor only", "Both"])
# plt2 = plt.bar([1,2,3], Deltier, bottom=Tratier, color='darkorange')
# plt3 = plt.bar([1,2,3], Duptier, bottom=Tratier+Deltier, color='darkorchid')
# plt4 = plt.bar([1,2,3], Invtier, bottom=Tratier+Deltier+Duptier, color='cornflowerblue')
# plt5 = plt.bar([1,2,3], Instier, bottom=Tratier+Deltier+Duptier+Invtier, color='saddlebrown')
# labels = ['TRA','DEL', 'DUP', 'INV', "INS"]
# ax.set_ylabel("Count")
# plt.legend(labels)
# plt.tight_layout()
# fig.savefig('countOccurrence_ifRestrictingLen.png')
# plt.close(fig)
#
# '''horizontal stacked barplot of the relative proportion of each sv type in sample type'''
# fig, ax = plt.subplots()
# Deltier = np.array((N_rel_Del, T_rel_Del, Both_rel_Del))
# Duptier = np.array((N_rel_Dup, T_rel_Dup, Both_rel_Dup))
# Invtier = np.array((N_rel_Inv, T_rel_Inv, Both_rel_Inv))
# Tratier = np.array((N_rel_Tra, T_rel_Tra, Both_rel_Tra))
# Instier = np.array((N_rel_Ins, T_rel_Ins, Both_rel_Ins))
# plt1 = plt.barh([1,2,3], Tratier, left=[0,0,0], color='darkturquoise', tick_label = ["Normal only", "Tumor only", "Both"])
# plt2 = plt.barh([1,2,3], Deltier, left=Tratier, color='darkorange')
# plt3 = plt.barh([1,2,3], Duptier, left=Tratier+Deltier, color='darkorchid')
# plt4 = plt.barh([1,2,3], Invtier, left=Tratier+Deltier+Duptier, color='cornflowerblue')
# plt5 = plt.barh([1,2,3], Instier, left=Tratier+Deltier+Duptier+Invtier, color='saddlebrown')
#
# labels = ['TRA','DEL', 'DUP', 'INV', "INS"]
# plt.legend(labels)
# plt.tight_layout()
# fig.savefig('relativeOccurrence_ifRestrictingLen.png')
# plt.close(fig)

'''histograms of the lengths split by sample type and sv type'''
fig, axes = plt.subplots(4,3, figsize=(30,30))
font = {'family' : 'arial'}
matplotlib.rc('font', **font)


cols = ['Normal Only', 'Tumor Only', 'Both']
rows = ['DEL', 'DUP', 'INV', 'INS']

pad = 5 # in points

for ax, col in zip(axes[0], cols):
    ax.annotate(col, xy=(0.5, 1), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size='50', ha='center', va='baseline')

for ax, row in zip(axes[:,0], rows):
    ax.annotate(row, xy=(0, 0.5), xytext=(-ax.yaxis.labelpad - pad, 0),
                xycoords=ax.yaxis.label, textcoords='offset points',
                size='50', ha='right', va='center')


axes[0,0].hist(N_type_length['DEL'], bins = 50, color = 'darkorange')
axes[0,0].set_yscale('log')
#axes[0,0].set_xscale('log')
axes[1,0].hist(N_type_length['DUP'], bins = 50, color = 'darkorchid')
axes[1,0].set_yscale('log')
#axes[1,0].set_xscale('log')
axes[2,0].hist(N_type_length['INV'], bins = 50, color = 'cornflowerblue')
axes[2,0].set_yscale('log')
#axes[2,0].set_xscale('log')
axes[3,0].hist(N_type_length['INS'], bins = 50, color = 'saddlebrown')
axes[3,0].set_yscale('log')
#axes[3,0].set_xscale('log')

axes[0,1].hist(T_type_length['DEL'], bins = 50, color = 'darkorange')
axes[0,1].set_yscale('log')
#axes[0,1].set_xscale('log')
axes[1,1].hist(T_type_length['DUP'], bins = 50, color = 'darkorchid')
axes[1,1].set_yscale('log')
axes[1,1].set_xticks([0, 10000, 20000, 30000, 40000, 50000])
axes[1,1].set_xticklabels(['0', '10000', '20000', '30000', '40000', '50000'])
#axes[1,1].set_xscale('log')
axes[2,1].hist(T_type_length['INV'], bins = 50, color = 'cornflowerblue')
axes[2,1].set_yscale('log')
#axes[2,1].set_xscale('log')
axes[3,1].hist(T_type_length['INS'], bins = 50, color = 'saddlebrown')
axes[3,1].set_yscale('log')
#axes[3,1].set_xscale('log')

axes[0,2].hist(Both_type_length['DEL'], bins = 50, color = 'darkorange')
axes[0,2].set_yscale('log')
#axes[0,2].set_xscale('log')
axes[1,2].hist(Both_type_length['DUP'], bins = 50, color = 'darkorchid')
axes[1,2].set_yscale('log')
#axes[1,2].set_xscale('log')
axes[2,2].hist(Both_type_length['INV'], bins = 50, color = 'cornflowerblue')
axes[2,2].set_yscale('log')
#axes[2,2].set_xscale('log')
axes[3,2].hist(Both_type_length['INS'], bins = 50, color = 'saddlebrown')
axes[3,2].set_yscale('log')
#axes[3,2].set_xscale('log')

for i in range(len(rows)):
    for j in range(len(cols)):
        axes[i,j].set_xlabel('Length (bp)', fontsize=30)
        axes[i,j].set_ylabel('Log(count)', fontsize=30)
        for tick in axes[i,j].xaxis.get_major_ticks():
            tick.label.set_fontsize(25)
        for tick in axes[i,j].yaxis.get_major_ticks():
            tick.label.set_fontsize(30)

fig.savefig("SVlengths_ylog_ifRestrictingLen_betterLabels.png")
plt.tight_layout()
plt.close(fig)

quit()
'''histograms of the lengths split by sample type and sv type, cut off lengths'''
fig, axes = plt.subplots(4,3, figsize=(20,20),sharex='col')
plt.setp(axes.flat, ylabel='Count')
cols = ['SVs in Normal Only', 'SVs in Tumor Only', 'SVs in Both']
rows = ['DEL', 'DUP', 'INV', 'INS']

pad = 5 # in points


for ax, col in zip(axes[0], cols):
    ax.annotate(col, xy=(0.5, 1), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size='xx-large', ha='center', va='baseline')


for ax, row in zip(axes[:,0], rows):
    ax.annotate(row, xy=(0, 0.5), xytext=(-ax.yaxis.labelpad - pad, 0),
                xycoords=ax.yaxis.label, textcoords='offset points',
                size='x-large', ha='right', va='center')


axes[0,0].hist(N_type_length['DEL'], range=(0,50000), bins = 50, color = 'darkorange')
axes[1,0].hist(N_type_length['DUP'], range=(0,50000), bins = 50, color = 'darkorchid')
axes[2,0].hist(N_type_length['INV'], range=(0,50000), bins = 50, color = 'cornflowerblue')
axes[3,0].hist(N_type_length['INS'], range=(0,50000), bins = 50, color = 'saddlebrown')
axes[3,0].set_xlabel('Length (bp)')

axes[0,1].hist(T_type_length['DEL'], range=(0,50000), bins = 50, color = 'darkorange')
axes[1,1].hist(T_type_length['DUP'], range=(0,50000), bins = 50, color = 'darkorchid')
axes[2,1].hist(T_type_length['INV'], range=(0,50000), bins = 50, color = 'cornflowerblue')
axes[3,1].hist(T_type_length['INS'], range=(0,50000), bins = 50, color = 'saddlebrown')
axes[3,1].set_xlabel('Length (bp)')

axes[0,2].hist(Both_type_length['DEL'], range=(0,50000), bins = 50, color = 'darkorange')
axes[1,2].hist(Both_type_length['DUP'], range=(0,50000), bins = 50, color = 'darkorchid')
axes[2,2].hist(Both_type_length['INV'], range=(0,50000), bins = 50, color = 'cornflowerblue')
axes[3,2].hist(Both_type_length['INS'], range=(0,50000), bins = 50, color = 'saddlebrown')
axes[3,2].set_xlabel('Length (bp)')

fig.savefig("SVlengths_limitLength_ifRestrictingLen.png")
plt.close(fig)

'''histograms of the lengths split by sample type and sv type, cut off lengths, and shared axes'''
fig, axes = plt.subplots(4,3, figsize=(20,20), sharey='row', sharex='col')
cols = ['SVs in Normal Only', 'SVs in Tumor Only', 'SVs in Both']
rows = ['DEL', 'DUP', 'INV', 'INS']

pad = 5 # in points


for ax, col in zip(axes[0], cols):
    ax.annotate(col, xy=(0.5, 1), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size='xx-large', ha='center', va='baseline')


for ax, row in zip(axes[:,0], rows):
    ax.annotate(row, xy=(0, 0.5), xytext=(-ax.yaxis.labelpad - pad, 0),
                xycoords=ax.yaxis.label, textcoords='offset points',
                size='x-large', ha='right', va='center')


axes[0,0].hist(N_type_length['DEL'], range=(0,50000), bins = 50, color = 'darkorange')
axes[0,0].set_ylabel('Count')
axes[1,0].hist(N_type_length['DUP'], range=(0,50000), bins = 50, color = 'darkorchid')
axes[1,0].set_ylabel('Count')
axes[2,0].hist(N_type_length['INV'], range=(0,50000), bins = 50, color = 'cornflowerblue')
axes[2,0].set_ylabel('Count')
axes[3,0].hist(N_type_length['INS'], range=(0,50000), bins = 50, color = 'saddlebrown')
axes[3,0].set_ylabel('Count')
axes[3,0].set_xlabel('Length (bp)')

axes[0,1].hist(T_type_length['DEL'], range=(0,50000), bins = 50, color = 'darkorange')
axes[1,1].hist(T_type_length['DUP'], range=(0,50000), bins = 50, color = 'darkorchid')
axes[2,1].hist(T_type_length['INV'], range=(0,50000), bins = 50, color = 'cornflowerblue')
axes[3,1].hist(T_type_length['INS'], range=(0,50000), bins = 50, color = 'saddlebrown')
axes[3,1].set_xlabel('Length (bp)')

axes[0,2].hist(Both_type_length['DEL'], range=(0,50000), bins = 50, color = 'darkorange')
axes[1,2].hist(Both_type_length['DUP'], range=(0,50000), bins = 50, color = 'darkorchid')
axes[2,2].hist(Both_type_length['INV'], range=(0,50000), bins = 50, color = 'cornflowerblue')
axes[3,2].hist(Both_type_length['INS'], range=(0,50000), bins = 50, color = 'saddlebrown')
axes[3,2].set_xlabel('Length (bp)')

fig.savefig("SVlengths_limitLength_shareAxes_ifRestrictingLen.png")
plt.close(fig)

'''boxplots of the lengths split by sample type and sv type'''
fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(15,10), sharey=True)
names = ['DEL', 'DUP', 'INV', 'INS']
vals_N = [N_type_length['DEL'], N_type_length['DUP'], N_type_length['INV'], N_type_length['INS']]
vals_T = [T_type_length['DEL'], T_type_length['DUP'], T_type_length['INV'], T_type_length['INS']]
vals_Both = [Both_type_length['DEL'], Both_type_length['DUP'], Both_type_length['INV'], Both_type_length['INS']]


#N
xs_N = []
for i, value in enumerate([N_num_Del, N_num_Dup, N_num_Inv, N_num_Ins]):
    xs_N.append(np.random.normal(i+1, 0.04, value))


#T
xs_T = []
for i, value in enumerate([T_num_Del, T_num_Dup, T_num_Inv, T_num_Ins]):
    xs_T.append(np.random.normal(i+1, 0.04, value))

#Both
xs_Both = []
for i, value in enumerate([Both_num_Del, Both_num_Dup, Both_num_Inv, Both_num_Ins]):
    xs_Both.append(np.random.normal(i+1, 0.04, value))

ymin = 0
ymax = 50000
ax1.set_title("SVs in Normal Only", size='xx-large')
ax1.set_ylabel("Length (bp)", size='x-large')
ax1.boxplot(vals_N, labels=names, boxprops = {'linewidth':3}, medianprops = {'linewidth':3})
ax1.set_ylim([ymin, ymax])
for xs, values in zip(xs_N, vals_N):
    ax1.scatter(xs, values, alpha=0.2)
ax2.set_title("SVs in Tumor Only", size='xx-large')
ax2.boxplot(vals_T, labels=names, boxprops = {'linewidth':3}, medianprops = {'linewidth':3})
ax2.set_ylim([ymin, ymax])
ax2.set_xlabel("SV type", size='x-large')
for xs, values in zip(xs_T, vals_T):
    ax2.scatter(xs, values, alpha=0.2)
ax3.set_title("SVs in Both", size='xx-large')
ax3.boxplot(vals_Both, labels=names, boxprops = {'linewidth':3}, medianprops = {'linewidth':3})
ax3.set_ylim([ymin, ymax])
for xs, values in zip(xs_Both, vals_Both):
    ax3.scatter(xs, values, alpha=0.2)

fig.savefig("boxplot_withScatter_ifRestrictingLen.png")
plt.close(fig)

'''histograms of the lengths split only by sample type'''
fig, (ax1, ax2, ax3) = plt.subplots(1,3)
ax1.set_title("SVs in Normal Only")
ax1.hist(N_lengths, bins=50, color='dodgerblue')
ax1.set_yscale('log')
ax2.set_title("SVs in Tumor Only")
ax2.hist(T_lengths, bins=50, color='olive')
ax2.set_yscale('log')
ax3.set_title("SVs in Both")
ax3.hist(Both_lengths, bins=50, color='saddlebrown')
ax3.set_yscale('log')
ax2.set_xlabel("Length (bp)")
plt.tight_layout()
fig.savefig("lengthbySample_ifRestrictingLen.png")
plt.close(fig)

'''boxplots of the lengths split only by sample type'''
fig, ax = plt.subplots()
names = ["Normal Only", "Tumor Only", "Both"]
ymin = 0
ymax = 50000

xs = []
for i, value in enumerate([len(N_lengths), len(T_lengths), len(Both_lengths)]):
    xs.append(np.random.normal(i+1, 0.04, value))

ax.set_ylim([ymin, ymax])
plt.boxplot([N_lengths, T_lengths, Both_lengths], labels=names, boxprops = {'linewidth':3}, medianprops = {'linewidth':3})
for xs, values in zip(xs, [N_lengths, T_lengths, Both_lengths]):
    plt.scatter(xs, values, alpha=0.2)
fig.savefig("boxplot_lengthsBySample_withScatter_ifRestrictingLen.png")
plt.close(fig)
