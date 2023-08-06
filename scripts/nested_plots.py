import argparse
import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import pysam
from collections import defaultdict
import mappy as mp
import matplotlib.font_manager as font_manager
import matplotlib
import os.path

def get_diff_list(list1, list2):
    ret = []
    for i in range(len(list1)):
        ret.append(list2[i]-list1[i])
    return ret


parser = argparse.ArgumentParser( description='Generate plots positional ambiguity before and after ')
parser.add_argument('--input-tsv', required=True)
parser.add_argument('--output-prefix', required=True)
parser.add_argument('--output-dir', required=True)
args = parser.parse_args()

alu_dups = {}
sva_dups = {}
combined_dups = {}
with open(args.input_tsv, 'r') as in_tsv:
    for line in in_tsv:
        row = line.strip().split('\t')
        if "Duplication" in row[1]:
            dv = float(row[1].split('_')[1])
            total_reads = int(row[4])
            if total_reads == 0:
                continue
            if dv not in combined_dups:
                combined_dups[dv] = defaultdict(list)
            combined_dups[dv][1000].append(1- int(row[5])/total_reads)
            combined_dups[dv][500].append(1- int(row[6])/total_reads)
            combined_dups[dv][250].append(1- int(row[7])/total_reads)
            combined_dups[dv][100].append(1- int(row[8])/total_reads)
            combined_dups[dv][50].append(1- int(row[9])/total_reads)
            if "Alu" in row[3]:
                if dv not in alu_dups:
                    alu_dups[dv] = defaultdict(list)
                alu_dups[dv][1000].append(1- int(row[5])/total_reads)
                alu_dups[dv][500].append(1- int(row[6])/total_reads)
                alu_dups[dv][250].append(1- int(row[7])/total_reads)
                alu_dups[dv][100].append(1- int(row[8])/total_reads)
                alu_dups[dv][50].append(1- int(row[9])/total_reads)
            elif "SVA" in row[3]:
                if dv not in sva_dups:
                    sva_dups[dv] = defaultdict(list)
                sva_dups[dv][1000].append(1- int(row[5])/total_reads)
                sva_dups[dv][500].append(1- int(row[6])/total_reads)
                sva_dups[dv][250].append(1- int(row[7])/total_reads)
                sva_dups[dv][100].append(1- int(row[8])/total_reads)
                sva_dups[dv][50].append(1- int(row[9])/total_reads)


# Compute the mean and stdev of each measure
data = defaultdict(list)
count = 1
for dv in combined_dups:
    data['X'].append(count)
    data['XString'].append(str(count -1))
    count += 1

data['X'] = np.array(data['X'])
wd = 0.4
    

for dv in combined_dups:
    data["combined_250"].append(combined_dups[dv][250])
    data["alu_250"].append(alu_dups[dv][250])
    data["sva_250"].append(sva_dups[dv][250])
    data["combined_50"].append(combined_dups[dv][50])
    data["alu_50"].append(alu_dups[dv][50])
    data["sva_50"].append(sva_dups[dv][50])


plt.boxplot(data["combined_250"], positions=data['X'], widths=wd)
plt.xticks(data['X'], data['XString'], fontsize=10)
plt.yticks(fontsize=15)
plt.title('Combined SVA & Alu Positional Ambiguity within +/- 250bp vs Sequence Divergence', fontsize=20)
plt.xlabel('Sequence Divergence %', fontsize=17)
plt.ylabel('Fraction of Reads with Positional Ambiguity +/-250bp', fontsize=17)
fig = plt.gcf()
fig.set_size_inches(24.5, 10.5)
fig.savefig(args.output_dir+"/"+args.output_prefix+"_combined_250.png")
plt.clf()


plt.boxplot(data["combined_50"], positions=data['X'], widths=wd)
plt.xticks(data['X'], data['XString'], fontsize=10)
plt.yticks(fontsize=15)
plt.title('Combined SVA & Alu Positional Ambiguity within +/- 50bp vs Sequence Divergence', fontsize=20)
plt.xlabel('Sequence Divergence %', fontsize=17)
plt.ylabel('Fraction of Reads with Positional Ambiguity +/-50bp', fontsize=17)
fig = plt.gcf()
fig.set_size_inches(24.5, 10.5)
fig.savefig(args.output_dir+"/"+args.output_prefix+"_combined_50.png")
plt.clf()


plt.boxplot(data["alu_250"], positions=data['X'], widths=wd)
plt.xticks(data['X'], data['XString'], fontsize=10)
plt.yticks(fontsize=15)
plt.title('Alu Positional Ambiguity within +/- 250bp vs Sequence Divergence', fontsize=20)
plt.xlabel('Sequence Divergence %', fontsize=17)
plt.ylabel('Fraction of Reads with Positional Ambiguity +/-250bp', fontsize=17)
fig = plt.gcf()
fig.set_size_inches(24.5, 10.5)
fig.savefig(args.output_dir+"/"+args.output_prefix+"_alu_250.png")
plt.clf()


plt.boxplot(data["alu_50"], positions=data['X'], widths=wd)
plt.xticks(data['X'], data['XString'], fontsize=14)
plt.yticks(fontsize=20)
plt.title('Misaligned Alu Insertions within +/- 50bp vs Sequence Divergence', fontsize=28)
plt.xlabel('Sequence Divergence %', fontsize=23)
plt.ylabel('Fraction of Reads with Positional Ambiguity +/-50bp', fontsize=23)
fig = plt.gcf()
fig.set_size_inches(28, 12.5)
fig.savefig(args.output_dir+"/"+args.output_prefix+"_alu_50.png")
plt.clf()

plt.boxplot(data["sva_250"], positions=data['X'], widths=wd)
plt.xticks(data['X'], data['XString'], fontsize=10)
plt.yticks(fontsize=15)
plt.title('SVA Positional Ambiguity within +/- 250bp vs Sequence Divergence', fontsize=20)
plt.xlabel('Sequence Divergence %', fontsize=17)
plt.ylabel('Fraction of Reads with Positional Ambiguity +/-250bp', fontsize=17)
fig = plt.gcf()
fig.set_size_inches(24.5, 10.5)
fig.savefig(args.output_dir+"/"+args.output_prefix+"_sva_250.png")
plt.clf()


plt.boxplot(data["sva_50"], positions=data['X'], widths=wd)
plt.xticks(data['X'], data['XString'], fontsize=10)
plt.yticks(fontsize=15)
plt.title('SVA Positional Ambiguity within +/- 50bp vs Sequence Divergence', fontsize=20)
plt.xlabel('Sequence Divergence %', fontsize=17)
plt.ylabel('Fraction of Reads with Positional Ambiguity +/-50bp', fontsize=17)
fig = plt.gcf()
fig.set_size_inches(24.5, 10.5)
fig.savefig(args.output_dir+"/"+args.output_prefix+"_sva_50.png")
plt.clf()


plt.subplot(2, 2, 1)
plt.boxplot(data["alu_250"], positions=data['X'], widths=wd)
plt.xticks(data['X'], data['XString'], fontsize=15)
plt.yticks(fontsize=25)
plt.ylim(top=1.01)
plt.title('Alu Positional Ambiguity within +/- 250bp vs Sequence Divergence', fontsize=25)
plt.xlabel('Sequence Divergence %', fontsize=20)
plt.ylabel('Fraction of Reads with Positional Ambiguity +/-250bp', fontsize=20)

plt.subplot(2, 2, 2)
plt.boxplot(data["alu_50"], positions=data['X'], widths=wd)
plt.xticks(data['X'], data['XString'], fontsize=15)
plt.yticks(fontsize=25)
plt.ylim(top=1.01)
plt.title('Alu Positional Ambiguity within +/- 50bp vs Sequence Divergence', fontsize=25)
plt.xlabel('Sequence Divergence %', fontsize=20)
plt.ylabel('Fraction of Reads with Positional Ambiguity +/-50bp', fontsize=20)

plt.subplot(2, 2, 3)
plt.boxplot(data["sva_250"], positions=data['X'], widths=wd)
plt.xticks(data['X'], data['XString'], fontsize=15)
plt.yticks(fontsize=25)
plt.ylim(top=1.01)
plt.title('SVA Positional Ambiguity within +/- 250bp vs Sequence Divergence', fontsize=25)
plt.xlabel('Sequence Divergence %', fontsize=20)
plt.ylabel('Fraction of Reads with Positional Ambiguity +/-250bp', fontsize=20)

plt.subplot(2, 2, 4)
plt.boxplot(data["sva_50"], positions=data['X'], widths=wd)
plt.xticks(data['X'], data['XString'], fontsize=15)
plt.yticks(fontsize=25)
plt.ylim(top=1.01)
plt.title('SVA Positional Ambiguity within +/- 50bp vs Sequence Divergence', fontsize=25)
plt.xlabel('Sequence Divergence %', fontsize=20)
plt.ylabel('Fraction of Reads with Positional Ambiguity +/-50bp', fontsize=20)

plt.tight_layout()
fig = plt.gcf()
fig.set_size_inches(55, 35)
fig.savefig(args.output_dir+"/"+args.output_prefix+"_sva_alu_all.png")
plt.clf()
