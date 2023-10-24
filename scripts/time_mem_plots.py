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


parser = argparse.ArgumentParser( description='Generate plots of before and after realignment data')
parser.add_argument('--input-dir', required=True)
parser.add_argument('--output-prefix', required=True)
parser.add_argument('--output-dir', required=True)
parser.add_argument('--samples', required=True)
parser.add_argument('--fractions', required=True)
args = parser.parse_args()

# Read in the data for what we're interested in
tldr_before_mem = defaultdict(list) 
tldr_before_time = defaultdict(list)
xtea_before_mem = defaultdict(list)
xtea_before_time = defaultdict(list)
somrit_mem = defaultdict(list)
somrit_time = defaultdict(list)
somrit_ideal_mem = defaultdict(list)
somrit_ideal_time = defaultdict(list)
somrit_realign_only_mem = defaultdict(list)
somrit_realign_only_time = defaultdict(list)

fractions = args.fractions.split(',')
samples = args.samples.split(',')

for sample in samples:
    for f in fractions:
        frac = int(f)
        for rep in [1,2,3]:
            # check if the file exists
            if not os.path.isfile(args.input_dir+"/"+sample+"/Time_Mem_"+sample+"_"+f+"_"+str(rep)+".txt"):
                continue
            with open(args.input_dir+"/"+sample+"/Time_Mem_"+sample+"_"+f+"_"+str(rep)+".txt", 'r') as in_time_mem:
                for line in in_time_mem:
                    row = line.strip().split("\t")
                    xtea_before_mem[f].append(float(row[7]))
                    xtea_before_time[f].append(float(row[8]))
                    tldr_before_mem[f].append(float(row[9]))
                    tldr_before_time[f].append(float(row[10]))
                    somrit_extract_mem = float(row[11])
                    somrit_extract_time= float(row[12])
                    somrit_realign_mem = float(row[13])
                    somrit_realign_time= float(row[14])
                    somrit_realign_max_mem = float(row[15])
                    somrit_realign_max_time= float(row[16])
                    somrit_realign_avg_mem = float(row[17])
                    somrit_realign_avg_time= float(row[18])
                    somrit_merge_mem = float(row[19])
                    somrit_merge_time= float(row[20])
                    somrit_classify_mem = float(row[21])
                    somrit_classify_time= float(row[22])
                    somrit_filter_mem = float(row[23])
                    somrit_filter_time= float(row[24])
                    somrit_time[f].append(somrit_extract_time + somrit_realign_time + somrit_merge_time + somrit_classify_time + somrit_filter_time)
                    somrit_mem[f].append(max(somrit_extract_mem, somrit_realign_max_mem, somrit_merge_mem, somrit_classify_mem, somrit_filter_mem))
                    somrit_realign_only_mem[f].append(somrit_realign_max_mem)
                    somrit_realign_only_time[f].append(somrit_realign_time)
                    somrit_ideal_mem[f].append(max(somrit_extract_mem, somrit_realign_max_mem, somrit_merge_mem, somrit_classify_mem, somrit_filter_mem))
                    somrit_ideal_time[f].append(somrit_extract_time + somrit_realign_max_time + somrit_merge_time + somrit_classify_time + somrit_filter_time)

print("data")

# Compute the mean and stdev of each measure
data = defaultdict(list)
count = 0
for cov in fractions:
    data['X'].append(int(cov))
    data['XString'].append(str(cov))
    data["xtea_before_mem_mean"].append(np.mean(xtea_before_mem[cov]))
    data["xtea_before_time_mean"].append(np.mean(xtea_before_time[cov]))
    data["tldr_before_mem_mean"].append(np.mean(tldr_before_mem[cov]))
    data["tldr_before_time_mean"].append(np.mean(tldr_before_time[cov]))
    data["somrit_mem_mean"].append(np.mean(somrit_mem[cov]))
    data["somrit_time_mean"].append(np.mean(somrit_time[cov]))
    data["somrit_realign_only_mem_mean"].append(np.mean(somrit_realign_only_mem[cov]))
    data["somrit_realign_only_time_mean"].append(np.mean(somrit_realign_only_time[cov]))
    data["somrit_ideal_mem_mean"].append(np.mean(somrit_ideal_mem[cov]))
    data["somrit_ideal_time_mean"].append(np.mean(somrit_ideal_time[cov]))
    data["xtea_before_mem"].append(xtea_before_mem[cov])
    data["xtea_before_time"].append(xtea_before_time[cov])
    data["tldr_before_mem"].append(tldr_before_mem[cov])
    data["tldr_before_time"].append(tldr_before_time[cov])
    data["somrit_mem"].append(somrit_mem[cov])
    data["somrit_time"].append(somrit_time[cov])
    data["somrit_ideal_mem"].append(somrit_ideal_mem[cov])
    data["somrit_ideal_time"].append(somrit_ideal_time[cov])
    data["somrit_realign_only_mem"].append(somrit_realign_only_mem[cov])
    data["somrit_realign_only_time"].append(somrit_realign_only_time[cov])

data['X'] = np.array(data['X'])
# Plot Memory Graph
wd = 0.4
plt.bar(data['X']+2*wd,data["xtea_before_mem_mean"], color='b', alpha=0.4, width=wd, edgecolor='k', label="xTea-Long")
plt.bar(data['X']+3*wd,data["tldr_before_mem_mean"], color='g', alpha=0.4, width=wd, edgecolor='k', label="tldr")
plt.bar(data['X']+4*wd,data["somrit_mem_mean"], color='c', alpha=0.4, width=wd, edgecolor='k', label="somrit")
plt.boxplot(data["xtea_before_mem"], positions=data['X']+2*wd, widths=wd)
plt.boxplot(data["tldr_before_mem"], positions=data['X']+3*wd, widths=wd)
plt.boxplot(data["somrit_mem"], positions=data['X']+4*wd, widths=wd)

plt.xticks(data['X']+2.5*wd, data['XString'], fontsize=15)
plt.yticks(fontsize=15)
plt.title('Maximum Memory Usage per Program', fontsize=20)
plt.xlabel('Coverage', fontsize=17)
plt.ylabel('Memory Usage (MiB)', fontsize=17)
plt.legend(loc='upper left', fontsize=15)
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
fig.savefig(args.output_dir+"/"+args.output_prefix+"_Memory.png")
plt.clf()


plt.bar(data['X']+2*wd,data["xtea_before_time_mean"], color='b', alpha=0.4, width=wd, edgecolor='k', label="xTea-Long")
plt.bar(data['X']+3*wd,data["tldr_before_time_mean"], color='g', alpha=0.4, width=wd, edgecolor='k', label="tldr")
plt.bar(data['X']+4*wd,data["somrit_time_mean"], color='c', alpha=0.4, width=wd, edgecolor='k', label="somrit total")
plt.bar(data['X']+5*wd,data["somrit_ideal_time_mean"], alpha=0.4, color='m', width=wd, edgecolor='k', label="somrit ideal")
plt.boxplot(data["xtea_before_time"], positions=data['X']+2*wd, widths=wd)
plt.boxplot(data["tldr_before_time"], positions=data['X']+3*wd, widths=wd)
plt.boxplot(data["somrit_time"], positions=data['X']+4*wd, widths=wd)
plt.boxplot(data["somrit_ideal_time"], positions=data['X']+5*wd, widths=wd)

plt.xticks(data['X']+2.5*wd, data['XString'], fontsize=15)
plt.yticks(fontsize=15)
plt.title('Total Time per Program', fontsize=20)
plt.xlabel('Coverage', fontsize=17)
plt.ylabel('Time (s)', fontsize=17)
plt.legend(loc='upper left', fontsize=15)
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
fig.savefig(args.output_dir+"/"+args.output_prefix+"_Time.png")
plt.clf()
print("Time")

