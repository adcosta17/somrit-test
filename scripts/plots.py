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
parser.add_argument('--somrit-dir', required=True)
parser.add_argument('--output-prefix', required=True)
parser.add_argument('--output-dir', required=True)
parser.add_argument('--samples', required=True)
parser.add_argument('--fractions', required=True)
parser.add_argument('--coverage', required=True)
args = parser.parse_args()

# Read in the data for what we're interested in
cute_sv_before_mem = defaultdict(list)
cute_sv_before_time = defaultdict(list)
sniffles_before_mem = defaultdict(list)
sniffles_before_time = defaultdict(list)
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

i = 0
fractions = args.fractions.split(':')
sample_coverage = {}
coverages = args.coverage.split(',')
sample_fractions = {}
all_coverage = []
for sample in args.samples.split(','):
    sample_fractions[sample] = fractions[i].split(',')
    sample_coverage[sample] = int(coverages[i])
    if i == 0:
        # poplate the all_coverage list
        for f in sample_fractions[sample]:
            all_coverage.append(round(float(f)*int(coverages[i])))
    i += 1

for sample in sample_fractions:
    for f in sample_fractions[sample]:
        frac = float(f)
        for rep in [1,2,3]:
            # check if the file exists
            if not os.path.isfile(args.input_dir+"/"+sample+"/Time_Mem_"+sample+"_"+f+"_"+str(rep)+".txt"):
                continue
            with open(args.input_dir+"/"+sample+"/Time_Mem_"+sample+"_"+f+"_"+str(rep)+".txt", 'r') as in_time_mem:
                for line in in_time_mem:
                    row = line.strip().split("\t")
                    cute_sv_before_mem[round(frac*sample_coverage[sample])].append(float(row[3]))
                    cute_sv_before_time[round(frac*sample_coverage[sample])].append(float(row[4]))
                    sniffles_before_mem[round(frac*sample_coverage[sample])].append(float(row[5]))
                    sniffles_before_time[round(frac*sample_coverage[sample])].append(float(row[6]))
                    xtea_before_mem[round(frac*sample_coverage[sample])].append(float(row[7]))
                    xtea_before_time[round(frac*sample_coverage[sample])].append(float(row[8]))
                    tldr_before_mem[round(frac*sample_coverage[sample])].append(float(row[9]))
                    tldr_before_time[round(frac*sample_coverage[sample])].append(float(row[10]))
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
                    somrit_time[round(frac*sample_coverage[sample])].append(somrit_extract_time + somrit_realign_time + somrit_merge_time + somrit_classify_time + somrit_filter_time)
                    somrit_mem[round(frac*sample_coverage[sample])].append(max(somrit_extract_mem, somrit_realign_max_mem, somrit_merge_mem, somrit_classify_mem, somrit_filter_mem))
                    somrit_realign_only_mem[round(frac*sample_coverage[sample])].append(somrit_realign_max_mem)
                    somrit_realign_only_time[round(frac*sample_coverage[sample])].append(somrit_realign_time)
                    somrit_ideal_mem[round(frac*sample_coverage[sample])].append(max(somrit_extract_mem, somrit_realign_max_mem, somrit_merge_mem, somrit_classify_mem, somrit_filter_mem))
                    somrit_ideal_time[round(frac*sample_coverage[sample])].append(somrit_extract_time + somrit_realign_max_time + somrit_merge_time + somrit_classify_time + somrit_filter_time)

print("data")

# Compute the mean and stdev of each measure
data = defaultdict(list)
count = 0
for cov in all_coverage:
    data['X'].append(cov)
    data['XString'].append(str(cov))
    data["cute_sv_before_mem_mean"].append(np.mean(cute_sv_before_mem[cov]))
    data["cute_sv_before_time_mean"].append(np.mean(cute_sv_before_time[cov]))
    data["sniffles_before_mem_mean"].append(np.mean(sniffles_before_mem[cov]))
    data["sniffles_before_time_mean"].append(np.mean(sniffles_before_time[cov]))
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
    data["cute_sv_before_mem"].append(cute_sv_before_mem[cov])
    data["cute_sv_before_time"].append(cute_sv_before_time[cov])
    data["sniffles_before_mem"].append(sniffles_before_mem[cov])
    data["sniffles_before_time"].append(sniffles_before_time[cov])
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
plt.bar(data['X'],data["cute_sv_before_mem_mean"], color='r', alpha=0.4, width=wd, edgecolor='k', label="CuteSV")
plt.bar(data['X']+wd,data["sniffles_before_mem_mean"], color='y', alpha=0.4, width=wd, edgecolor='k', label="Sniffles")
plt.bar(data['X']+2*wd,data["xtea_before_mem_mean"], color='b', alpha=0.4, width=wd, edgecolor='k', label="xTea")
plt.bar(data['X']+3*wd,data["tldr_before_mem_mean"], color='g', alpha=0.4, width=wd, edgecolor='k', label="tldr")
plt.bar(data['X']+4*wd,data["somrit_mem_mean"], color='c', alpha=0.4, width=wd, edgecolor='k', label="somrit")
plt.boxplot(data["cute_sv_before_mem"], positions=data['X'], widths=wd)
plt.boxplot(data["sniffles_before_mem"], positions=data['X']+wd, widths=wd)
plt.boxplot(data["xtea_before_mem"], positions=data['X']+2*wd, widths=wd)
plt.boxplot(data["tldr_before_mem"], positions=data['X']+3*wd, widths=wd)
plt.boxplot(data["somrit_mem"], positions=data['X']+4*wd, widths=wd)

plt.xticks(data['X']+2.5*wd, data['XString'], fontsize=15)
plt.yticks(fontsize=15)
plt.title('Max Memory Usage per Tool', fontsize=20)
plt.xlabel('Coverage', fontsize=17)
plt.ylabel('Memory Usage (MiB)', fontsize=17)
plt.legend(loc='upper left', fontsize=15)
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
fig.savefig(args.output_dir+"/"+args.output_prefix+"_Memory.png")
plt.clf()


plt.bar(data['X'],data["cute_sv_before_time_mean"], color='r', alpha=0.4, width=wd, edgecolor='k', label="CuteSV")
plt.bar(data['X']+wd,data["sniffles_before_time_mean"], color='y', alpha=0.4, width=wd, edgecolor='k', label="Sniffles")
plt.bar(data['X']+2*wd,data["xtea_before_time_mean"], color='b', alpha=0.4, width=wd, edgecolor='k', label="xTea")
plt.bar(data['X']+3*wd,data["tldr_before_time_mean"], color='g', alpha=0.4, width=wd, edgecolor='k', label="tldr")
plt.bar(data['X']+4*wd,data["somrit_time_mean"], color='c', alpha=0.4, width=wd, edgecolor='k', label="somrit total")
plt.bar(data['X']+5*wd,data["somrit_ideal_time_mean"], alpha=0.4, color='m', width=wd, edgecolor='k', label="somrit ideal")
plt.boxplot(data["cute_sv_before_time"], positions=data['X'], widths=wd)
plt.boxplot(data["sniffles_before_time"], positions=data['X']+wd, widths=wd)
plt.boxplot(data["xtea_before_time"], positions=data['X']+2*wd, widths=wd)
plt.boxplot(data["tldr_before_time"], positions=data['X']+3*wd, widths=wd)
plt.boxplot(data["somrit_time"], positions=data['X']+4*wd, widths=wd)
plt.boxplot(data["somrit_ideal_time"], positions=data['X']+5*wd, widths=wd)

plt.xticks(data['X']+2.5*wd, data['XString'], fontsize=15)
plt.yticks(fontsize=15)
plt.title('Max Time per Tool', fontsize=20)
plt.xlabel('Coverage', fontsize=17)
plt.ylabel('Time (s)', fontsize=17)
plt.legend(loc='upper left', fontsize=15)
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
fig.savefig(args.output_dir+"/"+args.output_prefix+"_Time.png")
plt.clf()
print("Time")

# Read and Plot FP Translocations before/after
cute_sv_before = defaultdict(list)
cute_sv_after = defaultdict(list)
sniffles_before = defaultdict(list)
sniffles_after = defaultdict(list)
min_sniffles = 100
max_sniffes = 0
min_cutesv = 100
max_cutesv = 0
for sample in sample_fractions:
    for f in sample_fractions[sample]:
        frac = float(f)
        for rep in [1,2,3]:
            # check if the file exists
            if not os.path.isfile(args.input_dir+"/"+sample.replace('0','')+"/"+sample.replace('0','')+"_Summary_translocations_"+f+"_"+str(rep)+".txt"):
                continue
            with open(args.input_dir+"/"+sample.replace('0','')+"/"+sample.replace('0','')+"_Summary_translocations_"+f+"_"+str(rep)+".txt", 'r') as in_fp_tra:
                for line in in_fp_tra:
                    row = line.strip().split("\t")
                    cute_sv_before[round(frac*sample_coverage[sample])].append(int(row[4]))
                    cute_sv_after[round(frac*sample_coverage[sample])].append(int(row[5]))
                    sniffles_before[round(frac*sample_coverage[sample])].append(int(row[2]))
                    sniffles_after[round(frac*sample_coverage[sample])].append(int(row[3]))
                    if (int(row[2])-int(row[3]))/int(row[2]) > max_sniffes:
                        max_sniffes = (int(row[2])-int(row[3]))/int(row[2])
                    if (int(row[2])-int(row[3]))/int(row[2]) < min_sniffles:
                        min_sniffles = (int(row[2])-int(row[3]))/int(row[2])
                    if (int(row[4])-int(row[5]))/int(row[4]) > max_cutesv:
                        max_cutesv = (int(row[4])-int(row[5]))/int(row[4])
                    if (int(row[4])-int(row[5]))/int(row[4]) < min_cutesv:
                        min_cutesv = (int(row[4])-int(row[5]))/int(row[4])

print(min_sniffles)
print(max_sniffes)
print(min_cutesv)
print(max_cutesv)

for cov in all_coverage:
    data["cute_sv_before_tra_mean"].append(np.mean(cute_sv_before[cov]))
    data["cute_sv_after_tra_mean"].append(np.mean(cute_sv_after[cov]))
    data["sniffles_before_tra_mean"].append(np.mean(sniffles_before[cov]))
    data["sniffles_after_tra_mean"].append(np.mean(sniffles_after[cov]))
    data["cute_sv_before_tra"].append(cute_sv_before[cov])
    data["cute_sv_after_tra"].append(cute_sv_after[cov])
    data["sniffles_before_tra"].append(sniffles_before[cov])
    data["sniffles_after_tra"].append(sniffles_after[cov])

wd = 0.4
plt.bar(data['X'],data["cute_sv_before_tra_mean"], color='r', alpha=0.4, width=wd, edgecolor='k', label="CuteSV Before")
plt.bar(data['X']+wd,data["cute_sv_after_tra_mean"], color='y', alpha=0.4, width=wd, edgecolor='k', label="CuteSV After")
plt.bar(data['X']+2*wd+0.2,data["sniffles_before_tra_mean"], color='b', alpha=0.4, width=wd, edgecolor='k', label="Sniffles2 Before")
plt.bar(data['X']+3*wd+0.2,data["sniffles_after_tra_mean"], color='g', alpha=0.4, width=wd, edgecolor='k', label="Sniffles2 After")
plt.boxplot(data["cute_sv_before_tra"], positions=data['X'], widths=wd)
plt.boxplot(data["cute_sv_after_tra"], positions=data['X']+wd, widths=wd)
plt.boxplot(data["sniffles_before_tra"], positions=data['X']+2*wd+0.2, widths=wd)
plt.boxplot(data["sniffles_after_tra"], positions=data['X']+3*wd+0.2, widths=wd)

plt.xticks(data['X']+1.5*wd, data['XString'], fontsize=15)
plt.yticks(fontsize=15)
plt.title('False Positive Translocations Before and After Realignment', fontsize=20)
plt.xlabel('Coverage', fontsize=17)
plt.ylabel('FP Translocations', fontsize=17)
plt.legend(loc='upper left', fontsize=15)
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
fig.savefig(args.output_dir+"/"+args.output_prefix+"_FP_TRA.png")
plt.clf()


# Read and Plot Sniffles/CuteSV Before after Precision & Recall
xtea_precision = defaultdict(list)
tldr_precision = defaultdict(list)
somrit_precision = defaultdict(list)
xtea_recall = defaultdict(list)
tldr_recall = defaultdict(list)
somrit_recall = defaultdict(list)
for sample in sample_fractions:
    #print(sample)
    for f in sample_fractions[sample]:
        #print(f)
        frac = float(f)
        for rep in [1,2,3]:
            #print(rep)
            # check if the file exists
            if not os.path.isfile(args.input_dir+"/"+sample+"/Somrit_xtea_tldr_before_after_"+sample.replace('0','')+"_"+f+"_"+str(rep)+".txt"):
                print(args.input_dir+"/"+sample+"/Somrit_xtea_tldr_before_after_"+sample.replace('0','')+"_"+f+"_"+str(rep)+".txt")
                continue
            with open(args.input_dir+"/"+sample+"/Somrit_xtea_tldr_before_after_"+sample.replace('0','')+"_"+f+"_"+str(rep)+".txt", 'r') as in_fp_tra:
                for line in in_fp_tra:
                    row = line.strip().split("\t")
                    if float(row[6]) > 0.1:
                        xtea_precision[round(frac*sample_coverage[sample])].append(float(row[5]))
                        xtea_recall[round(frac*sample_coverage[sample])].append(float(row[6]))
                    tldr_precision[round(frac*sample_coverage[sample])].append(float(row[10]))
                    tldr_recall[round(frac*sample_coverage[sample])].append(float(row[11]))
            if not os.path.isfile(args.somrit_dir+"/"+sample+"/Somrit_xtea_tldr_before_after_"+sample.replace('0','')+"_"+f+"_"+str(rep)+".txt"):
                continue
            with open(args.somrit_dir+"/"+sample+"/Somrit_xtea_tldr_before_after_"+sample.replace('0','')+"_"+f+"_"+str(rep)+".txt", 'r') as in_fp_tra:
                for line in in_fp_tra:
                    somrit_precision[round(frac*sample_coverage[sample])].append(float(row[15]))
                    somrit_recall[round(frac*sample_coverage[sample])].append(float(row[16]))



for cov in all_coverage:
    data["xtea_precision_mean"].append(np.mean(xtea_precision[cov]))
    data["tldr_precision_mean"].append(np.mean(tldr_precision[cov]))
    data["somrit_precision_mean"].append(np.mean(somrit_precision[cov]))
    data["xtea_precision"].append(xtea_precision[cov])
    data["tldr_precision"].append(tldr_precision[cov])
    data["somrit_precision"].append(somrit_precision[cov])
    data["xtea_recall_mean"].append(np.mean(xtea_recall[cov]))
    data["tldr_recall_mean"].append(np.mean(tldr_recall[cov]))
    data["somrit_recall_mean"].append(np.mean(somrit_recall[cov]))
    data["xtea_recall"].append(xtea_recall[cov])
    data["tldr_recall"].append(tldr_recall[cov])
    data["somrit_recall"].append(somrit_recall[cov])


wd = 0.4
plt.bar(data['X'],data["xtea_precision_mean"], color='r', width=wd, edgecolor='k', label="Xtea Precision", alpha=0.4)
plt.bar(data['X']+2*wd,data["tldr_precision_mean"], color='b', width=wd, edgecolor='k', label="Tldr Precision", alpha=0.4)
plt.bar(data['X']+4*wd,data["somrit_precision_mean"], color='c', width=wd, edgecolor='k', label="Somrit Precision", alpha=0.4)
plt.bar(data['X']+wd,data["xtea_recall_mean"], color='y', width=wd, edgecolor='k', label="Xtea Recall", alpha=0.4)
plt.bar(data['X']+3*wd,data["tldr_recall_mean"], color='g', width=wd, edgecolor='k', label="Tldr Recall", alpha=0.4)
plt.bar(data['X']+5*wd,data["somrit_recall_mean"], color='m', width=wd, edgecolor='k', label="Somrit Recall", alpha=0.4)
plt.boxplot(data["xtea_precision"], positions=data['X'], widths=wd)
plt.boxplot(data["tldr_precision"], positions=data['X']+2*wd, widths=wd)
plt.boxplot(data["somrit_precision"], positions=data['X']+4*wd, widths=wd)
plt.boxplot(data["xtea_recall"], positions=data['X']+wd, widths=wd)
plt.boxplot(data["tldr_recall"], positions=data['X']+3*wd, widths=wd)
plt.boxplot(data["somrit_recall"], positions=data['X']+5*wd, widths=wd)

plt.xticks(data['X']+2.5*wd, data['XString'], fontsize=15)
plt.yticks(fontsize=15)
plt.title('Somrit, Xtea and Tldr Precision and Recall', fontsize=20)
plt.xlabel('Coverage', fontsize=17)
plt.ylabel('Value', fontsize=17)
plt.legend(loc='upper left', fontsize=15)
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
fig.savefig(args.output_dir+"/"+args.output_prefix+"_Somrit_Xtea_Tldr_PR.png")
plt.clf()

wd = 0.8
plt.bar(data['X'],data["xtea_precision_mean"], color='r', width=wd, edgecolor='k', label="Xtea Precision", alpha=0.4)
plt.bar(data['X']+wd,data["tldr_precision_mean"], color='b', width=wd, edgecolor='k', label="Tldr Precision", alpha=0.4)
plt.bar(data['X']+2*wd,data["somrit_precision_mean"], color='c', width=wd, edgecolor='k', label="Somrit Precision", alpha=0.4)
plt.boxplot(data["xtea_precision"], positions=data['X'], widths=wd)
plt.boxplot(data["tldr_precision"], positions=data['X']+wd, widths=wd)
plt.boxplot(data["somrit_precision"], positions=data['X']+2*wd, widths=wd)

plt.xticks(data['X']+1*wd, data['XString'], fontsize=15)
plt.yticks(fontsize=15)
plt.title('Somrit, Xtea and Tldr Precision', fontsize=20)
plt.xlabel('Coverage', fontsize=17)
plt.ylabel('Precision', fontsize=17)
plt.legend(loc='upper left', fontsize=15)
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
fig.savefig(args.output_dir+"/"+args.output_prefix+"_Somrit_Xtea_Tldr_Precision.png")
plt.clf()

wd = 0.8
plt.bar(data['X'],data["xtea_recall_mean"], color='y', width=wd, edgecolor='k', label="Xtea Recall", alpha=0.4)
plt.bar(data['X']+wd,data["tldr_recall_mean"], color='g', width=wd, edgecolor='k', label="Tldr Recall", alpha=0.4)
plt.bar(data['X']+2*wd,data["somrit_recall_mean"], color='m', width=wd, edgecolor='k', label="Somrit Recall", alpha=0.4)
plt.boxplot(data["xtea_recall"], positions=data['X'], widths=wd)
plt.boxplot(data["tldr_recall"], positions=data['X']+wd, widths=wd)
plt.boxplot(data["somrit_recall"], positions=data['X']+2*wd, widths=wd)

plt.xticks(data['X']+1*wd, data['XString'], fontsize=15)
plt.yticks(fontsize=15)
plt.title('Somrit, Xtea and Tldr Recall', fontsize=20)
plt.xlabel('Coverage', fontsize=17)
plt.ylabel('Recall', fontsize=17)
plt.legend(loc='upper left', fontsize=15)
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
fig.savefig(args.output_dir+"/"+args.output_prefix+"_Somrit_Xtea_Tldr_Recall.png")
plt.clf()


plt.subplot(2, 1, 1)
plt.bar(data['X'],data["xtea_precision_mean"], color='r', width=wd, edgecolor='k', label="Xtea-Long", alpha=0.4)
plt.bar(data['X']+wd,data["tldr_precision_mean"], color='b', width=wd, edgecolor='k', label="Tldr", alpha=0.4)
plt.bar(data['X']+2*wd,data["somrit_precision_mean"], color='c', width=wd, edgecolor='k', label="Somrit", alpha=0.4)
plt.boxplot(data["xtea_precision"], positions=data['X'], widths=wd)
plt.boxplot(data["tldr_precision"], positions=data['X']+wd, widths=wd)
plt.boxplot(data["somrit_precision"], positions=data['X']+2*wd, widths=wd)
plt.xticks(data['X']+1*wd, data['XString'], fontsize=22)
plt.yticks(fontsize=22)
plt.ylim(top=0.95)
plt.title('Somrit, Xtea-Long and Tldr Precision', fontsize=27)
plt.xlabel('Coverage', fontsize=24)
plt.ylabel('Precision', fontsize=24)
plt.legend(loc='lower right', fontsize=19)

plt.subplot(2, 1, 2)
plt.bar(data['X'],data["xtea_recall_mean"], color='r', width=wd, edgecolor='k', label="Xtea-Long", alpha=0.4)
plt.bar(data['X']+wd,data["tldr_recall_mean"], color='b', width=wd, edgecolor='k', label="Tldr", alpha=0.4)
plt.bar(data['X']+2*wd,data["somrit_recall_mean"], color='g', width=wd, edgecolor='k', label="Somrit", alpha=0.4)
plt.boxplot(data["xtea_recall"], positions=data['X'], widths=wd)
plt.boxplot(data["tldr_recall"], positions=data['X']+wd, widths=wd)
plt.boxplot(data["somrit_recall"], positions=data['X']+2*wd, widths=wd)
plt.xticks(data['X']+1*wd, data['XString'], fontsize=22)
plt.yticks(fontsize=22)
plt.ylim(top=0.95)
plt.title('Somrit, Xtea-Long and Tldr Recall', fontsize=27)
plt.xlabel('Coverage', fontsize=24)
plt.ylabel('Recall', fontsize=24)
plt.legend(loc='lower right', fontsize=19)

plt.tight_layout()
fig = plt.gcf()
fig.set_size_inches(16.5, 22.5)
fig.savefig(args.output_dir+"/"+args.output_prefix+"_Somrit_Xtea_Tldr_Combined_Precision_Recall.png")
plt.clf()


plt.subplot(1, 2, 1)
plt.bar(data['X'],data["xtea_precision_mean"], color='r', width=wd, edgecolor='k', label="Xtea-Long", alpha=0.4)
plt.bar(data['X']+wd,data["tldr_precision_mean"], color='b', width=wd, edgecolor='k', label="Tldr", alpha=0.4)
plt.bar(data['X']+2*wd,data["somrit_precision_mean"], color='g', width=wd, edgecolor='k', label="Somrit", alpha=0.4)
plt.boxplot(data["xtea_precision"], positions=data['X'], widths=wd)
plt.boxplot(data["tldr_precision"], positions=data['X']+wd, widths=wd)
plt.boxplot(data["somrit_precision"], positions=data['X']+2*wd, widths=wd)
plt.xticks(data['X']+1*wd, data['XString'], fontsize=20)
plt.yticks(fontsize=20)
plt.ylim(top=0.95)
plt.title('Somrit, Xtea-Long and Tldr Precision', fontsize=25)
plt.xlabel('Coverage', fontsize=22)
plt.ylabel('Precision', fontsize=22)
plt.legend(loc='lower right', fontsize=17)

plt.subplot(1, 2, 2)
plt.bar(data['X'],data["xtea_recall_mean"], color='r', width=wd, edgecolor='k', label="Xtea-Long", alpha=0.4)
plt.bar(data['X']+wd,data["tldr_recall_mean"], color='b', width=wd, edgecolor='k', label="Tldr", alpha=0.4)
plt.bar(data['X']+2*wd,data["somrit_recall_mean"], color='g', width=wd, edgecolor='k', label="Somrit", alpha=0.4)
plt.boxplot(data["xtea_recall"], positions=data['X'], widths=wd)
plt.boxplot(data["tldr_recall"], positions=data['X']+wd, widths=wd)
plt.boxplot(data["somrit_recall"], positions=data['X']+2*wd, widths=wd)
plt.ylim(top=0.95)
plt.xticks(data['X']+1*wd, data['XString'], fontsize=20)
plt.yticks(fontsize=20)
plt.title('Somrit, Xtea-Long and Tldr Recall', fontsize=25)
plt.xlabel('Coverage', fontsize=22)
plt.ylabel('Recall', fontsize=22)
plt.legend(loc='lower right', fontsize=17)

plt.tight_layout()
fig = plt.gcf()
fig.set_size_inches(35, 12.5)
fig.savefig(args.output_dir+"/"+args.output_prefix+"_Somrit_Xtea_Tldr_Combined_Precision_Recall_Beside.png")
plt.clf()
