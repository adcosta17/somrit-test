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
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition, mark_inset)

def get_diff_list(list1, list2):
    ret = []
    for i in range(len(list1)):
        ret.append(list2[i]-list1[i])
    return ret


parser = argparse.ArgumentParser( description='Generate plots of before and after realignment data')
parser.add_argument('--input-dir', required=True)
parser.add_argument('--output-prefix', required=True)
parser.add_argument('--output-dir', required=True)
parser.add_argument('--samples', required=True)
args = parser.parse_args()

data = defaultdict(list)
# Read and Plot Sniffles/CuteSV Before after Precision & Recall
xtea_recall = defaultdict(list)
sniffles_recall = defaultdict(list)
tldr_recall = defaultdict(list)
somrit_recall = defaultdict(list)
tldr_read_recall = defaultdict(list)
somrit_read_recall = defaultdict(list)
somrit_precision = defaultdict(list)
tldr_precision = defaultdict(list)
for sample in args.samples.split(','):
    for rep in [1,2,3,4,5,6,7,8,9,10,11,12]:
        if not os.path.isfile(args.input_dir+"/"+sample+"/simulation_results_"+sample.replace('0','')+".spike_in."+str(rep)+".txt"):
            print("Missing: " + args.input_dir+"/"+sample+"/simulation_results_"+sample.replace('0','')+".spike_in."+str(rep)+".txt")
            continue
        with open(args.input_dir+"/"+sample+"/simulation_results_"+sample.replace('0','')+".spike_in."+str(rep)+".txt", 'r') as in_tsv:
            for line in in_tsv:
                row = line.strip().split("\t")
                somrit_recall[sample].append(float(row[3]))
                somrit_read_recall[sample].append(float(row[6]))
                tldr_recall[sample].append(float(row[9]))
                tldr_read_recall[sample].append(float(row[12]))
                xtea_recall[sample].append(float(row[15]))
                sniffles_recall[sample].append(float(row[18]))
                somrit_precision[sample].append(float(row[23]))
                print(row[23] + "\t" + row[25])
                tldr_precision[sample].append(float(row[25]))
                break

i = 1
wd = 0.5
count = 0
for sample in args.samples.split(','):
    data["X"].append(float(i))
    data["X1"].append(float(i)+wd)
    data["X2"].append(float(i)+2*wd)
    data["X3"].append(float(i)+3*wd)
    data["X1.5"].append(float(i)+1.5*wd)
    data["XString"].append(sample)
    data["sniffles_mean"].append(np.mean(sniffles_recall[sample]))
    data["sniffles"].append(sniffles_recall[sample])
    data["somrit_mean"].append(np.mean(somrit_recall[sample]))
    data["somrit"].append(somrit_recall[sample])
    data["somrit_read_mean"].append(np.mean(somrit_read_recall[sample]))
    data["somrit_read"].append(somrit_read_recall[sample])
    data["xtea_mean"].append(np.mean(xtea_recall[sample]))
    data["xtea"].append(xtea_recall[sample])
    data["tldr_mean"].append(np.mean(tldr_recall[sample]))
    data["tldr"].append(tldr_recall[sample])
    data["tldr_read_mean"].append(np.mean(tldr_read_recall[sample]))
    data["tldr_read"].append(tldr_read_recall[sample])
    data["tldr_precision"].append(tldr_precision[sample])
    data["tldr_precision_mean"].append(np.mean(tldr_precision[sample]))
    data["somrit_precision"].append(somrit_precision[sample])
    data["somrit_precision_mean"].append(np.mean(somrit_precision[sample]))
    i += 2.5


fig, ax1 = plt.subplots()
wd = 0.5
ax1.bar(data['X'],data["somrit_mean"], color='c', width=wd, edgecolor='k', label="Somrit", alpha=0.4)
ax1.bar(data['X1'],data["tldr_mean"], color='b', width=wd, edgecolor='k', label="tldr", alpha=0.4)
#plt.bar(data['X3'],data["tldr_read_mean"], color='g', width=wd, edgecolor='k', label="tldr Read Recall", alpha=0.4)
ax1.bar(data['X2'],data["sniffles_mean"], color='m', width=wd, edgecolor='k', label="Sniffles2", alpha=0.4)
#ax1.bar(data['X3'],data["xtea_mean"], color='c', width=wd, edgecolor='k', label="xTea-Long Recall", alpha=0.4)
ax1.boxplot(data["somrit"], positions=data['X'], widths=wd)
ax1.boxplot(data["tldr"], positions=data['X1'], widths=wd)
#plt.boxplot(data["tldr_read"], positions=data['X3'], widths=wd)
#ax1.boxplot(data["xtea"], positions=data['X3'], widths=wd)
ax1.boxplot(data["sniffles"], positions=data['X2'], widths=wd)
ax1.legend(loc='lower right', fontsize=15)

ax1.set_xticks(data['X1'])
ax1.set_xticklabels(data['XString'], fontsize=18)
ax1.set_yticks([0.0,0.2,0.4,0.6,0.8,1.0])
ax1.set_yticklabels(["0.0","0.2","0.4","0.6","0.8","1.0"],fontsize=18)
ax1.tick_params(axis='both', which='minor', size=15)
ax1.set_title('Somrit, tldr and Sniffles2 Recall on Simulated Insertions', fontsize=22)
ax1.set_xlabel('Sample', fontsize=20)
ax1.set_ylabel('Recall', fontsize=20)
ax1.set_ylim([0,1.05])

ax2 = plt.axes([0,0,1,1])
# Manually set the position and relative size of the inset axes within ax1
ip = InsetPosition(ax1, [0.67,0.67,0.3,0.3])
ax2.set_axes_locator(ip)
ax2.bar(data['X'],data["somrit_read_mean"], color='y', width=wd, edgecolor='k', label="Somrit Read Recall", alpha=0.4)
ax2.boxplot(data["somrit_read"], positions=data['X'], widths=wd)
ax2.set_xticks(data['X'])
ax2.set_xticklabels(data['XString'])
ax2.legend(loc=0)
ax2.set_ylim([0,1.05])

fig.set_size_inches(18.5, 10.5)
fig.savefig(args.output_dir+"/"+args.output_prefix+"_recall.png")
plt.clf()


fig, ax1 = plt.subplots()
wd = 1
ax1.bar(data['X'],data["somrit_precision_mean"], color='c', width=wd, edgecolor='k', label="Somrit", alpha=0.4)
ax1.bar(data['X1'],data["tldr_precision_mean"], color='b', width=wd, edgecolor='k', label="tldr", alpha=0.4)
ax1.boxplot(data["somrit_precision"], positions=data['X'], widths=wd)
ax1.boxplot(data["tldr_precision"], positions=data['X1'], widths=wd)
ax1.legend(loc='lower right', fontsize=15)

ax1.set_xticks(data['X1'])
ax1.set_xticklabels(data['XString'], fontsize=18)
ax1.tick_params(axis='both', which='minor', size=15)
ax1.set_title('Somrit and tldr Precision on Simulated Insertions', fontsize=22)
ax1.set_xlabel('Sample', fontsize=20)
ax1.set_ylabel('Recall', fontsize=20)
ax1.set_ylim([0,1.05])

fig.set_size_inches(18.5, 10.5)
fig.savefig(args.output_dir+"/"+args.output_prefix+"_precision.png")
plt.clf()



