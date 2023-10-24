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
args = parser.parse_args()

fractions = args.fractions.split(',')

data = defaultdict(list)
count = 0
for cov in fractions:
    data['X'].append(int(cov))
    data['XString'].append(str(cov))

data['X'] = np.array(data['X'])

# Read and Plot Sniffles/CuteSV Before after Precision & Recall
xtea_precision = defaultdict(list)
tldr_precision = defaultdict(list)
somrit_precision = defaultdict(list)
xtea_recall = defaultdict(list)
tldr_recall = defaultdict(list)
somrit_recall = defaultdict(list)
for sample in args.samples.split(','):
    #print(sample)
    for f in fractions:
        #print(f)
        frac = int(f)
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
                        xtea_precision[f].append(float(row[5]))
                        xtea_recall[f].append(float(row[6]))
                    tldr_precision[f].append(float(row[10]))
                    tldr_recall[f].append(float(row[11]))
                    somrit_precision[f].append(float(row[15]))
                    somrit_recall[f].append(float(row[16]))



for cov in fractions:
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
plt.bar(data['X']+2*wd,data["somrit_precision_mean"], color='g', width=wd, edgecolor='k', label="Somrit", alpha=0.4)
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
