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
    print(len(list1))
    print(len(list2))
    for i in range(len(list1)):
        ret.append(list2[i]-list1[i])
    return ret


parser = argparse.ArgumentParser( description='Generate plots of before and after realignment data')
parser.add_argument('--input-dir', required=True)
parser.add_argument('--output-prefix', required=True)
parser.add_argument('--output-dir', required=True)
parser.add_argument('--samples', required=True)
parser.add_argument('--fractions', required=True)
parser.add_argument('--coverage', required=True)
args = parser.parse_args()

# Read in the data for what we're interested in

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

# Compute the mean and stdev of each measure
data = defaultdict(list)
count = 0
for cov in all_coverage:
    data['X'].append(cov)
    data['XString'].append(str(cov))

data['X'] = np.array(data['X'])

# Read and Plot Sniffles/CuteSV Before after Precision & Recall
sniffles_before_diff = defaultdict(list)
sniffles_after_diff = defaultdict(list)
cute_sv_before_diff = defaultdict(list)
cute_sv_after_diff = defaultdict(list)
sniffles_diff = defaultdict(list)
cute_sv_diff = defaultdict(list)
sniffles_before_reads = defaultdict(list)
sniffles_after_reads = defaultdict(list)
cute_sv_before_reads = defaultdict(list)
cute_sv_after_reads = defaultdict(list)
sniffles_reads = defaultdict(list)
cute_sv_reads = defaultdict(list)
for sample in sample_fractions:
    #print(sample)
    for f in sample_fractions[sample]:
        #print(f)
        frac = float(f)
        for rep in [1,2,3,4,5]:
            #print(rep)
            # check if the file exists
            if not os.path.isfile(args.input_dir+"/"+sample+"/Diff_Sniffles_CuteSV_before_after_"+sample.replace('0','')+"_"+f+"_"+str(rep)+".txt"):
                continue
            with open(args.input_dir+"/"+sample+"/Diff_Sniffles_CuteSV_before_after_"+sample.replace('0','')+"_"+f+"_"+str(rep)+".txt", 'r') as in_fp_tra:
                for line in in_fp_tra:
                    row = line.strip().split("\t")
                    sniffles_before_diff[round(frac*sample_coverage[sample])].extend(list(map(float,row[4].split(','))))
                    sniffles_after_diff[round(frac*sample_coverage[sample])].extend(list(map(float,row[5].split(','))))
                    cute_sv_before_diff[round(frac*sample_coverage[sample])].extend(list(map(float,row[8].split(','))))
                    cute_sv_after_diff[round(frac*sample_coverage[sample])].extend(list(map(float,row[9].split(','))))
                    #sniffles_diff[round(frac*sample_coverage[sample])].extend(get_diff_list(list(map(float,row[4].split(','))), list(map(float,row[5].split(',')))))
                    #cute_sv_diff[round(frac*sample_coverage[sample])].extend(get_diff_list(list(map(float,row[8].split(','))), list(map(float,row[9].split(',')))))
                    sniffles_before_reads[round(frac*sample_coverage[sample])].extend(list(map(int,row[2].split(','))))
                    sniffles_after_reads[round(frac*sample_coverage[sample])].extend(list(map(int,row[3].split(','))))
                    cute_sv_before_reads[round(frac*sample_coverage[sample])].extend(list(map(int,row[6].split(','))))
                    cute_sv_after_reads[round(frac*sample_coverage[sample])].extend(list(map(int,row[7].split(','))))
                    sniffles_reads[round(frac*sample_coverage[sample])].extend(get_diff_list(list(map(int,row[2].split(','))), list(map(int,row[3].split(',')))))
                    cute_sv_reads[round(frac*sample_coverage[sample])].extend(get_diff_list(list(map(int,row[6].split(','))), list(map(int,row[7].split(',')))))
                    


for cov in all_coverage:
    data["cute_sv_before_diff_mean"].append(np.median(cute_sv_before_diff[cov]))
    data["cute_sv_after_diff_mean"].append(np.median(cute_sv_after_diff[cov]))
    data["sniffles_before_diff_mean"].append(np.median(sniffles_before_diff[cov]))
    data["sniffles_after_diff_mean"].append(np.median(sniffles_after_diff[cov]))
    data["cute_sv_before_diff"].append(cute_sv_before_diff[cov])
    data["cute_sv_after_diff"].append(cute_sv_after_diff[cov])
    data["sniffles_before_diff"].append(sniffles_before_diff[cov])
    data["sniffles_after_diff"].append(sniffles_after_diff[cov])
    #data["cute_sv_diff"].append(cute_sv_diff[cov])
    #data["sniffles_diff"].append(sniffles_diff[cov])
    #data["cute_sv_diff_mean"].append(np.mean(cute_sv_diff[cov]))
    #data["sniffles_diff_mean"].append(np.mean(sniffles_diff[cov]))
    data["cute_sv_before_reads_mean"].append(np.median(cute_sv_before_reads[cov]))
    data["cute_sv_after_reads_mean"].append(np.median(cute_sv_after_reads[cov]))
    data["sniffles_before_reads_mean"].append(np.median(sniffles_before_reads[cov]))
    data["sniffles_after_reads_mean"].append(np.median(sniffles_after_reads[cov]))
    data["cute_sv_before_reads"].append(cute_sv_before_reads[cov])
    data["cute_sv_after_reads"].append(cute_sv_after_reads[cov])
    data["sniffles_before_reads"].append(sniffles_before_reads[cov])
    data["sniffles_after_reads"].append(sniffles_after_reads[cov])
    data["cute_sv_reads"].append(cute_sv_reads[cov])
    data["sniffles_reads"].append(sniffles_reads[cov])
    data["cute_sv_reads_mean"].append(np.median(cute_sv_reads[cov]))
    data["sniffles_reads_mean"].append(np.median(sniffles_reads[cov]))


wd = 0.5
plt.bar(data['X'],data["sniffles_before_diff_mean"], color='r', width=wd, edgecolor='k', label="Sniffles2 Before Diff", alpha=0.4)
plt.bar(data['X']+wd,data["sniffles_after_diff_mean"], color='y', width=wd, edgecolor='k', label="Sniffles2 After Diff", alpha=0.4)
plt.boxplot(data["sniffles_before_diff"], positions=data['X'], widths=wd)
plt.boxplot(data["sniffles_after_diff"], positions=data['X']+wd, widths=wd)
plt.bar(data['X']+2*wd,data["cute_sv_before_diff_mean"], color='b', width=wd, edgecolor='k', label="CuteSV Before Diff", alpha=0.4)
plt.bar(data['X']+3*wd,data["cute_sv_after_diff_mean"], color='g', width=wd, edgecolor='k', label="CuteSV After Diff", alpha=0.4)
plt.boxplot(data["cute_sv_before_diff"], positions=data['X']+2*wd, widths=wd)
plt.boxplot(data["cute_sv_after_diff"], positions=data['X']+3*wd, widths=wd)

plt.xticks(data['X']+1.5*wd, data['XString'], fontsize=15)
plt.yticks(fontsize=15)
plt.title('Sniffles2 & CuteSV Positional Difference Before and After Realignment', fontsize=20)
plt.xlabel('Coverage', fontsize=17)
plt.ylabel('Distance (bp)', fontsize=17)
plt.legend(loc='upper left', fontsize=15)
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
fig.savefig(args.output_dir+"/"+args.output_prefix+"_Sniffles_CuteSV_Diff.png")
plt.clf()

wd = 0.5
#plt.bar(data['X'],data["sniffles_diff_mean"], color='r', width=wd, edgecolor='k', label="Sniffles Diff", alpha=0.4)
#plt.bar(data['X']+wd,data["cute_sv_diff_mean"], color='y', width=wd, edgecolor='k', label="CuteSV Diff", alpha=0.4)
#plt.boxplot(data["sniffles_diff"], positions=data['X'], widths=wd)
#plt.boxplot(data["cute_sv_diff"], positions=data['X']+wd, widths=wd)

#plt.xticks(data['X']+0.5*wd, data['XString'], fontsize=15)
#plt.yticks(fontsize=15)
#plt.title('CuteSV & Sniffles Positional Difference After Realignment', fontsize=20)
#plt.xlabel('Coverage', fontsize=17)
#plt.ylabel('Distance (bp)', fontsize=17)
#plt.legend(loc='upper left', fontsize=15)
#fig = plt.gcf()
#fig.set_size_inches(18.5, 10.5)
#fig.savefig(args.output_dir+"/"+args.output_prefix+"_Overall_Diff.png")
#plt.clf()


wd = 0.5
plt.bar(data['X'],data["sniffles_before_reads_mean"], color='r', width=wd, edgecolor='k', label="Sniffles2 Before Read Support", alpha=0.4)
plt.bar(data['X']+wd,data["sniffles_after_reads_mean"], color='y', width=wd, edgecolor='k', label="Sniffles2 After Read Support", alpha=0.4)
plt.boxplot(data["sniffles_before_reads"], positions=data['X'], widths=wd)
plt.boxplot(data["sniffles_after_reads"], positions=data['X']+wd, widths=wd)
plt.bar(data['X']+2*wd,data["cute_sv_before_reads_mean"], color='b', width=wd, edgecolor='k', label="CuteSV Before Read Support", alpha=0.4)
plt.bar(data['X']+3*wd,data["cute_sv_after_reads_mean"], color='g', width=wd, edgecolor='k', label="CuteSV After Read Support", alpha=0.4)
plt.boxplot(data["cute_sv_before_reads"], positions=data['X']+2*wd, widths=wd)
plt.boxplot(data["cute_sv_after_reads"], positions=data['X']+3*wd, widths=wd)

plt.xticks(data['X']+1.5*wd, data['XString'], fontsize=15)
plt.yticks(fontsize=15)
plt.title('Sniffles2 & CuteSV Read Support Before and After Realignment', fontsize=20)
plt.xlabel('Coverage', fontsize=17)
plt.ylabel('Read Support', fontsize=17)
plt.legend(loc='upper left', fontsize=15)
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
fig.savefig(args.output_dir+"/"+args.output_prefix+"_Sniffles_CuteSV_Read_Support.png")
plt.clf()

wd = 0.5
plt.bar(data['X'],data["sniffles_reads_mean"], color='r', width=wd, edgecolor='k', label="Sniffles Read Support", alpha=0.4)
plt.bar(data['X']+wd,data["cute_sv_reads_mean"], color='y', width=wd, edgecolor='k', label="CuteSV Read Support", alpha=0.4)
#plt.boxplot(data["sniffles_reads"], positions=data['X'], widths=wd)
#plt.boxplot(data["cute_sv_reads"], positions=data['X']+wd, widths=wd)

plt.xticks(data['X']+0.5*wd, data['XString'], fontsize=15)
plt.yticks(fontsize=15)
plt.title('CuteSV & Sniffles Positional Difference After Realignment', fontsize=20)
plt.xlabel('Coverage', fontsize=17)
plt.ylabel('Read Support', fontsize=17)
plt.legend(loc='upper left', fontsize=15)
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
fig.savefig(args.output_dir+"/"+args.output_prefix+"_Overall_Read_Support.png")
plt.clf()


# Read and Plot Sniffles/CuteSV Before after Precision & Recall
sniffles_before_diff = defaultdict(list)
sniffles_after_diff = defaultdict(list)
cute_sv_before_diff = defaultdict(list)
cute_sv_after_diff = defaultdict(list)
sniffles_diff = defaultdict(list)
cute_sv_diff = defaultdict(list)
sniffles_before_reads = defaultdict(list)
sniffles_after_reads = defaultdict(list)
cute_sv_before_reads = defaultdict(list)
cute_sv_after_reads = defaultdict(list)
sniffles_reads = defaultdict(list)
cute_sv_reads = defaultdict(list)
for sample in sample_fractions:
    #print(sample)
    for f in sample_fractions[sample]:
        #print(f)
        frac = float(f)
        for rep in [1,2,3,4,5]:
            #print(rep)
            # check if the file exists
            if not os.path.isfile(args.input_dir+"/"+sample+"/RT_Diff_Sniffles_CuteSV_before_after_"+sample.replace('0','')+"_"+f+"_"+str(rep)+".txt"):
                continue
            with open(args.input_dir+"/"+sample+"/RT_Diff_Sniffles_CuteSV_before_after_"+sample.replace('0','')+"_"+f+"_"+str(rep)+".txt", 'r') as in_fp_tra:
                for line in in_fp_tra:
                    row = line.strip().split("\t")
                    sniffles_before_diff[round(frac*sample_coverage[sample])].extend(list(map(float,row[4].split(','))))
                    sniffles_after_diff[round(frac*sample_coverage[sample])].extend(list(map(float,row[5].split(','))))
                    #cute_sv_before_diff[round(frac*sample_coverage[sample])].extend(list(map(float,row[8].split(','))))
                    #cute_sv_after_diff[round(frac*sample_coverage[sample])].extend(list(map(float,row[9].split(','))))
                    #sniffles_diff[round(frac*sample_coverage[sample])].extend(get_diff_list(list(map(float,row[4].split(','))), list(map(float,row[5].split(',')))))
                    #cute_sv_diff[round(frac*sample_coverage[sample])].extend(get_diff_list(list(map(float,row[8].split(','))), list(map(float,row[9].split(',')))))
                    sniffles_before_reads[round(frac*sample_coverage[sample])].extend(list(map(int,row[2].split(','))))
                    sniffles_after_reads[round(frac*sample_coverage[sample])].extend(list(map(int,row[3].split(','))))
                    cute_sv_before_reads[round(frac*sample_coverage[sample])].extend(list(map(int,row[6].split(','))))
                    cute_sv_after_reads[round(frac*sample_coverage[sample])].extend(list(map(int,row[7].split(','))))
                    sniffles_reads[round(frac*sample_coverage[sample])].extend(get_diff_list(list(map(int,row[2].split(','))), list(map(int,row[3].split(',')))))
                    cute_sv_reads[round(frac*sample_coverage[sample])].extend(get_diff_list(list(map(int,row[6].split(','))), list(map(int,row[7].split(',')))))

for cov in all_coverage:
    #data["large_cute_sv_diff"].append(cute_sv_diff[cov])
    #data["large_sniffles_diff"].append(sniffles_diff[cov])
    #data["large_cute_sv_diff_mean"].append(np.median(cute_sv_diff[cov]))
    #data["large_sniffles_diff_mean"].append(np.median(sniffles_diff[cov]))
    data["large_cute_sv_before_diff_mean"].append(np.median(cute_sv_before_diff[cov]))
    data["large_cute_sv_after_diff_mean"].append(np.median(cute_sv_after_diff[cov]))
    data["large_sniffles_before_diff_mean"].append(np.median(sniffles_before_diff[cov]))
    data["large_sniffles_after_diff_mean"].append(np.median(sniffles_after_diff[cov]))
    data["large_cute_sv_before_diff"].append(cute_sv_before_diff[cov])
    data["large_cute_sv_after_diff"].append(cute_sv_after_diff[cov])
    data["large_sniffles_before_diff"].append(sniffles_before_diff[cov])
    data["large_sniffles_after_diff"].append(sniffles_after_diff[cov])
    data["large_cute_sv_before_reads_mean"].append(np.median(cute_sv_before_reads[cov]))
    data["large_cute_sv_after_reads_mean"].append(np.median(cute_sv_after_reads[cov]))
    data["large_sniffles_before_reads_mean"].append(np.median(sniffles_before_reads[cov]))
    data["large_sniffles_after_reads_mean"].append(np.median(sniffles_after_reads[cov]))
    data["large_cute_sv_before_reads"].append(cute_sv_before_reads[cov])
    data["large_cute_sv_after_reads"].append(cute_sv_after_reads[cov])
    data["large_sniffles_before_reads"].append(sniffles_before_reads[cov])
    data["large_sniffles_after_reads"].append(sniffles_after_reads[cov])
    data["large_cute_sv_reads"].append(cute_sv_reads[cov])
    data["large_sniffles_reads"].append(sniffles_reads[cov])
    data["large_cute_sv_reads_mean"].append(np.median(cute_sv_reads[cov]))
    data["large_sniffles_reads_mean"].append(np.median(sniffles_reads[cov]))


wd = 0.5
plt.bar(data['X'],data["large_sniffles_before_diff_mean"], color='r', width=wd, edgecolor='k', label="Sniffles2 Before Diff", alpha=0.4)
plt.bar(data['X']+wd,data["large_sniffles_after_diff_mean"], color='y', width=wd, edgecolor='k', label="Sniffles2 After Diff", alpha=0.4)
plt.boxplot(data["large_sniffles_before_diff"], positions=data['X'], widths=wd)
plt.boxplot(data["large_sniffles_after_diff"], positions=data['X']+wd, widths=wd)
plt.bar(data['X']+2*wd,data["large_cute_sv_before_diff_mean"], color='b', width=wd, edgecolor='k', label="CuteSV Before Diff", alpha=0.4)
plt.bar(data['X']+3*wd,data["large_cute_sv_after_diff_mean"], color='g', width=wd, edgecolor='k', label="CuteSV After Diff", alpha=0.4)
plt.boxplot(data["large_cute_sv_before_diff"], positions=data['X']+2*wd, widths=wd)
plt.boxplot(data["large_cute_sv_after_diff"], positions=data['X']+3*wd, widths=wd)

plt.xticks(data['X']+1.5*wd, data['XString'], fontsize=15)
plt.yticks(fontsize=15)
plt.title('Sniffles2 & CuteSV Positional Difference Before and After Realignment', fontsize=20)
plt.xlabel('Coverage', fontsize=17)
plt.ylabel('Distance (bp)', fontsize=17)
plt.legend(loc='upper left', fontsize=15)
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
fig.savefig(args.output_dir+"/"+args.output_prefix+"_Sniffles_CuteSV_large_Diff.png")
plt.clf()

#wd = 0.5
#plt.bar(data['X'],data["large_sniffles_diff_mean"], color='r', width=wd, edgecolor='k', label="Sniffles Diff", alpha=0.4)
#plt.bar(data['X']+wd,data["large_cute_sv_diff_mean"], color='y', width=wd, edgecolor='k', label="CuteSV Diff", alpha=0.4)
#plt.boxplot(data["large_sniffles_diff"], positions=data['X'], widths=wd)
#plt.boxplot(data["large_cute_sv_diff"], positions=data['X']+wd, widths=wd)

#plt.xticks(data['X']+0.5*wd, data['XString'], fontsize=15)
#plt.yticks(fontsize=15)
#plt.title('CuteSV & Sniffles Positional Difference After Realignment', fontsize=20)
#plt.xlabel('Coverage', fontsize=17)
#plt.ylabel('Distance (bp)', fontsize=17)
#plt.legend(loc='upper left', fontsize=15)
#fig = plt.gcf()
#fig.set_size_inches(18.5, 10.5)
#fig.savefig(args.output_dir+"/"+args.output_prefix+"_Overall_large_Diff.png")
#plt.clf()


wd = 0.5
plt.bar(data['X'],data["large_sniffles_before_reads_mean"], color='r', width=wd, edgecolor='k', label="Sniffles2 Before Read Support", alpha=0.4)
plt.bar(data['X']+wd,data["large_sniffles_after_reads_mean"], color='y', width=wd, edgecolor='k', label="Sniffles2 After Read Support", alpha=0.4)
plt.boxplot(data["large_sniffles_before_reads"], positions=data['X'], widths=wd)
plt.boxplot(data["large_sniffles_after_reads"], positions=data['X']+wd, widths=wd)
plt.bar(data['X']+2*wd,data["large_cute_sv_before_reads_mean"], color='b', width=wd, edgecolor='k', label="CuteSV Before Read Support", alpha=0.4)
plt.bar(data['X']+3*wd,data["large_cute_sv_after_reads_mean"], color='g', width=wd, edgecolor='k', label="CuteSV After Read Support", alpha=0.4)
plt.boxplot(data["large_cute_sv_before_reads"], positions=data['X']+2*wd, widths=wd)
plt.boxplot(data["large_cute_sv_after_reads"], positions=data['X']+3*wd, widths=wd)

plt.xticks(data['X']+1.5*wd, data['XString'], fontsize=15)
plt.yticks(fontsize=15)
plt.title('Sniffles2 & CuteSV Read Support Before and After Realignment', fontsize=20)
plt.xlabel('Coverage', fontsize=17)
plt.ylabel('Read Support', fontsize=17)
plt.legend(loc='upper left', fontsize=15)
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
fig.savefig(args.output_dir+"/"+args.output_prefix+"_Sniffles_CuteSV_Read_Support_large.png")
plt.clf()

wd = 0.5
plt.bar(data['X'],data["large_sniffles_reads_mean"], color='r', width=wd, edgecolor='k', label="Sniffles Read Support", alpha=0.4)
plt.bar(data['X']+wd,data["large_cute_sv_reads_mean"], color='y', width=wd, edgecolor='k', label="CuteSV Read Support", alpha=0.4)
#plt.boxplot(data["large_sniffles_reads"], positions=data['X'], widths=wd)
#plt.boxplot(data["large_cute_sv_reads"], positions=data['X']+wd, widths=wd)

plt.xticks(data['X']+0.5*wd, data['XString'], fontsize=15)
plt.yticks(fontsize=15)
plt.title('CuteSV & Sniffles Positional Difference After Realignment', fontsize=20)
plt.xlabel('Coverage', fontsize=17)
plt.ylabel('Read Support', fontsize=17)
plt.legend(loc='upper left', fontsize=15)
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
fig.savefig(args.output_dir+"/"+args.output_prefix+"_Overall_Read_Support_large.png")
plt.clf()