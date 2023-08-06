import pysam
import argparse
import sys
import csv
from collections import defaultdict
import gzip
import re
import json
from intervaltree import Interval, IntervalTree
from statistics import mean

def plot(name, x, y1, y2):
    plt.ylim(0, 1)
    plt.plot(x, y1, 'bo', linestyle = 'dotted')
    plt.plot(x, y2, 'r+', linestyle = 'dotted')
    fig = plt.gcf()
    fig.savefig(name)
    plt.clf()

def plot_bar(name, x, y1):
    plt.bar(x, height=y1*100)
    fig = plt.gcf()
    fig.savefig(name)
    plt.clf()

def get_truth_dict(file, window):
    truth_dict = defaultdict(IntervalTree)
    with open(file) as in_truth:
        for line in in_truth:
            row = line.strip().split('\t')
            if len(row) < 2:
                continue
            chrom = row[0]
            start = int(row[1])
            key = start
            truth_dict[chrom][(start-window):(start+window)] = key
    return truth_dict

def get_pb_truth_dict(file, window, chroms, min_size=None, max_size=None):
    base_tree = defaultdict(IntervalTree)
    truth_dict = defaultdict(IntervalTree)
    with open(file) as in_truth:
        for line in in_truth:
            row = line.strip().split('\t')
            if min_size is not None:
                if len(row[7]) < min_size:
                    continue
            if max_size is not None:
                if len(row[7]) > max_size:
                    continue
            chrom = row[0]
            if chrom not in chroms:
                continue
            start = int(row[1])-window
            end = int(row[2])+window
            key = start
            base_tree[chrom][start:end] = key
    #print("Pass2")
    #count = 0
    with open(file) as in_truth:
        for line in in_truth:
            row = line.strip().split('\t')
            #count += 1
            #if count % 10000 == 0:
                #print(count)
            if len(row[7]) < 100:
                continue
            chrom = row[0]
            if chrom not in chroms:
                continue
            start = int(row[1])
            end = int(row[2])
            key = start
            nearby = base_tree[chrom][start:end]
            nearby_truth = truth_dict[chrom][start:end]
            if len(nearby) > 2 and len(nearby_truth) == 0:
                truth_dict[chrom][(start-window):(end+window)] = key
    return truth_dict

def get_xtea_stats(file, truth_dict):
    seen = {}
    tp = 0
    fp = {}
    fn = 0
    tp_distance = 0
    for chrom in truth_dict:
        seen[chrom] = defaultdict(int)
    with open(file) as in_calls:
        for line in in_calls:
            row = line.strip().split('\t')
            chrom = row[0]
            start = int(row[1])
            # Check if the position is within window size of a truth call
            nearby = truth_dict[chrom][start:start+1]
            if len(nearby) > 0:
                # Have a tp
                for item in nearby:
                    #print(item.data)
                    if int(item.data) in seen[chrom]:
                        seen[chrom][int(item.data)] = min(abs(int(item.data)-start),seen[chrom][int(item.data)])
                    else:
                        seen[chrom][int(item.data)] = abs(int(item.data)-start)
            else:
                # Have a fp
                fp[chrom+":"+str(start)] = 1
        for chrom in truth_dict:
            for item in truth_dict[chrom]:
                if int(item.data) not in seen[chrom]:
                    fn += 1
        for chrom in seen:
            for pos in seen[chrom]:
                tp += 1
                tp_distance += seen[chrom][pos]
    return tp, fp, fn, tp_distance/tp


def get_tldr_stats(file, repeat, truth_dict):
    seen = {}
    tp = 0
    fp = 1
    fn = 0
    for chrom in truth_dict:
        seen[chrom] = defaultdict(int)
    count = 0
    with open(file) as in_calls:
        for line in in_calls:
            if count == 0:
                count += 1
                continue
            row = line.strip().split('\t')
            chrom = row[1]
            start = int(row[2])
            end = int(row[3])
            # Check if the position is within window size of a truth call
            if repeat in row[5]:
                nearby = truth_dict[chrom][start:end]
                if len(nearby) > 0:
                    # Have a tp
                    for item in nearby:
                        #print(item.data)
                        seen[chrom][int(item.data)] = 1
                else:
                    # Have a fp
                    fp += 1
        for chrom in truth_dict:
            for item in truth_dict[chrom]:
                if int(item.data) not in seen[chrom]:
                    fn += 1
        for chrom in seen:
            for pos in seen[chrom]:
                tp += 1
    return tp, fp, fn

def get_all(file, truth_dict):
    seen = {}
    tp = 0
    fp = 0
    fn = 0
    for chrom in truth_dict:
        seen[chrom] = defaultdict(int)
    with open(file) as in_calls:
        for line in in_calls:
            row = line.strip().split('\t')
            chrom = row[0]
            start = int(row[1])
            # Check if the position is within window size of a truth call
            nearby = truth_dict[chrom][start:start+1]
            if len(nearby) > 0:
                # Have a tp
                for item in nearby:
                    #print(item.data)
                    seen[chrom][int(item.data)] = 1
            else:
                # Have a fp
                fp += 1
        for chrom in truth_dict:
            for item in truth_dict[chrom]:
                if int(item.data) not in seen[chrom]:
                    fn += 1
        for chrom in seen:
            for pos in seen[chrom]:
                tp += 1
    return seen


def get_diff(before, after):
    counts = {}
    for chrom in before:
        if chrom not in counts:
            counts[chrom] = defaultdict(int)
        for pos in before[chrom]:
            counts[chrom][pos] += 1
    for chrom in after:
        if chrom not in counts:
            counts[chrom] = defaultdict(int)
        for pos in after[chrom]:
            counts[chrom][pos] += 2
    before_only = {}
    after_only = {}
    for chrom in counts:
        for pos in counts[chrom]:
            if counts[chrom][pos] == 3:
                continue
            if counts[chrom][pos] == 1:
                before_only[str(chrom)+":"+str(pos)] = 1
            if counts[chrom][pos] == 2:
                after_only[str(chrom)+":"+str(pos)] = 1
    print("before_only")
    for item in before_only:
        print(item)
    print("after_only")
    for item in after_only:
        print(item)


def merge_truth(truth_lines, truth_sines, truth_svas):
    truth_dict = defaultdict(IntervalTree)
    for chrom in truth_lines:
        for item in truth_lines[chrom]:
            start = item.begin
            end = item.end
            data = item.data
            truth_dict[chrom][start:end] = data
    for chrom in truth_sines:
        for item in truth_sines[chrom]:
            start = item.begin
            end = item.end
            data = item.data
            truth_dict[chrom][start:end] = data
    for chrom in truth_svas:
        for item in truth_svas[chrom]:
            start = item.begin
            end = item.end
            data = item.data
            truth_dict[chrom][start:end] = data
    return truth_dict

parser = argparse.ArgumentParser( description='Get metrics for xtea and tldr before and after runs')
parser.add_argument('--truth-line', required=True)
parser.add_argument('--truth-sine', required=True)
parser.add_argument('--truth-sva', required=True)
parser.add_argument('--truth-pb', required=True)
parser.add_argument('--replicates', type=int, default=5)
parser.add_argument('--line-xtea', default="classified_results.txt.merged_LINE1.txt")
parser.add_argument('--sine-xtea', default="classified_results.txt.merged_ALU.txt")
parser.add_argument('--sva-xtea', default="classified_results.txt.merged_SVA.txt")
parser.add_argument('--tldr-suffix', default='.table.txt')
parser.add_argument('--window-size', type=int, default=1000)
parser.add_argument('--chrom-list', default="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY")
args = parser.parse_args()


# Read in truth data for each repeat type
truth_lines = get_truth_dict(args.truth_line, args.window_size)
truth_sines = get_truth_dict(args.truth_sine, args.window_size)
truth_svas = get_truth_dict(args.truth_sva, args.window_size)
truth_overall = merge_truth(truth_lines, truth_sines, truth_svas)
truth_pb = get_pb_truth_dict(args.truth_pb, args.window_size, args.chrom_list.split(','))

xtea_before_res= {}
xtea_before_res["ALL"] = {}
xtea_before_res["PB"] = {}
xtea_after_res= {}
xtea_after_res["ALL"] = {}
xtea_after_res["PB"] = {}
xtea_before_res["ALL"]["Precision"] = []
xtea_before_res["ALL"]["Recall"] = []
xtea_before_res["PB"]["Precision"] = []
xtea_before_res["PB"]["Recall"] = []
xtea_after_res["ALL"]["Precision"] = []
xtea_after_res["ALL"]["Recall"] = []
xtea_after_res["PB"]["Precision"] = []
xtea_after_res["PB"]["Recall"] = []

#print("Source\tFraction\tCoverage\tRepliacte\tRepeat\tBeforeTP\tBeforeFP\tBeforeFN\tBeforePrecision\tBeforeRecall\tAfterTP\tAfterFP\tAfterFN\tAfterPrecision\tAfterRecall")
for cov in ["0.05", "0.1", "0.15", "0.2", "0.25", "0.3", "0.35", "0.4", "0.45", "0.5", "0.55", "0.6", "0.65", "0.7", "0.75", "0.8", "0.85", "0.9", "0.95"]:
    print(cov)
    xtea_before_res["ALL"][cov] = defaultdict(list)
    xtea_before_res["PB"][cov] = defaultdict(list)
    xtea_after_res["ALL"][cov] = defaultdict(list)
    xtea_after_res["PB"][cov] = defaultdict(list)
    for rep in range(args.replicates):
        rep = rep + 1
        print(rep)
        # Xtea
        tp, fp, fn, tp_distance = get_xtea_stats("HG002_"+cov+"_"+str(rep)+"/"+args.line_xtea, truth_overall)
        xtea_before_res["ALL"][cov]["Precision"].append(tp/(tp+len(fp)))
        xtea_before_res["ALL"][cov]["Recall"].append(tp/(tp+fn))
        tp, fp, fn, tp_distance = get_xtea_stats("HG002_"+cov+"_"+str(rep)+"/"+args.line_xtea, truth_pb)
        xtea_before_res["PB"][cov]["Precision"].append(tp/(tp+len(fp)))
        xtea_before_res["PB"][cov]["Recall"].append(tp/(tp+fn))
        tp, fp, fn, tp_distance = get_xtea_stats("Realigned_HG002_"+cov+"_"+str(rep)+"/"+args.line_xtea, truth_overall)
        xtea_after_res["ALL"][cov]["Precision"].append(tp/(tp+len(fp)))
        xtea_after_res["ALL"][cov]["Recall"].append(tp/(tp+fn))
        tp, fp, fn, tp_distance = get_xtea_stats("Realigned_HG002_"+cov+"_"+str(rep)+"/"+args.line_xtea, truth_pb)
        xtea_after_res["PB"][cov]["Precision"].append(tp/(tp+len(fp)))
        xtea_after_res["PB"][cov]["Recall"].append(tp/(tp+fn))

        # tldr
        #before_stats["tldr_LINE"]["tp"], before_stats["tldr_LINE"]["fp"], before_stats["tldr_LINE"]["fn"] = get_tldr_stats("tldr_HG002_"+cov+"_"+str(rep)+"/"+"HG002_"+cov+"_"+str(rep)+args.tldr_suffix, "LINE", truth_lines)
        #before_stats["tldr_SINE"]["tp"], before_stats["tldr_SINE"]["fp"], before_stats["tldr_SINE"]["fn"] = get_tldr_stats("tldr_HG002_"+cov+"_"+str(rep)+"/"+"HG002_"+cov+"_"+str(rep)+args.tldr_suffix, "SINE", truth_sines)
        #before_stats["tldr_SVA"]["tp"], before_stats["tldr_SVA"]["fp"], before_stats["tldr_SVA"]["fn"] = get_tldr_stats("tldr_HG002_"+cov+"_"+str(rep)+"/"+"HG002_"+cov+"_"+str(rep)+args.tldr_suffix, "SVA", truth_svas)
        #after_stats["tldr_LINE"]["tp"], after_stats["tldr_LINE"]["fp"], after_stats["tldr_LINE"]["fn"] = get_tldr_stats("tldr_Realigned_HG002_"+cov+"_"+str(rep)+"/"+"Realigned_HG002_"+cov+"_"+str(rep)+args.tldr_suffix, "LINE",  truth_lines)
        #after_stats["tldr_SINE"]["tp"], after_stats["tldr_SINE"]["fp"], after_stats["tldr_SINE"]["fn"] = get_tldr_stats("tldr_Realigned_HG002_"+cov+"_"+str(rep)+"/"+"Realigned_HG002_"+cov+"_"+str(rep)+args.tldr_suffix, "SINE", truth_sines)
        #after_stats["tldr_SVA"]["tp"], after_stats["tldr_SVA"]["fp"], after_stats["tldr_SVA"]["fn"] = get_tldr_stats("tldr_Realigned_HG002_"+cov+"_"+str(rep)+"/"+"Realigned_HG002_"+cov+"_"+str(rep)+args.tldr_suffix, "SVA", truth_svas)
    # Average out the replicates for Precision and Recall
    xtea_before_res["ALL"]["Precision"].append(mean(xtea_before_res["ALL"][cov]["Precision"]))
    xtea_before_res["ALL"]["Recall"].append(mean(xtea_before_res["ALL"][cov]["Recall"]))
    xtea_before_res["PB"]["Precision"].append(mean(xtea_before_res["PB"][cov]["Precision"]))
    xtea_before_res["PB"]["Recall"].append(mean(xtea_before_res["PB"][cov]["Recall"]))
    xtea_after_res["ALL"]["Precision"].append(mean(xtea_after_res["ALL"][cov]["Precision"]))
    xtea_after_res["ALL"]["Recall"].append(mean(xtea_after_res["ALL"][cov]["Recall"]))
    xtea_after_res["PB"]["Precision"].append(mean(xtea_after_res["PB"][cov]["Precision"]))
    xtea_after_res["PB"]["Recall"].append(mean(xtea_after_res["PB"][cov]["Recall"]))

# Plot Before and after for each repeat family
import matplotlib.pyplot as plt
import numpy as np
print("Precision")
x = np.array(["3", "6", "9", "12", "15", "18", "21", "24", "27", "30", "33", "36", "39", "42", "45", "48", "51", "54", "57"])
plot('/u/adcosta/transfer/xtea_tldr_plots/ALL_Precision.png', x, np.array(xtea_before_res["ALL"]["Precision"]), np.array(xtea_after_res["ALL"]["Precision"]))
plot('/u/adcosta/transfer/xtea_tldr_plots/PB_Precision.png', x, np.array(xtea_before_res["PB"]["Precision"]), np.array(xtea_after_res["PB"]["Precision"]))
print("Recall")
plot('/u/adcosta/transfer/xtea_tldr_plots/ALL_Recall.png', x, np.array(xtea_before_res["ALL"]["Recall"]), np.array(xtea_after_res["ALL"]["Recall"]))
plot('/u/adcosta/transfer/xtea_tldr_plots/PB_Recall.png', x, np.array(xtea_before_res["PB"]["Recall"]), np.array(xtea_after_res["PB"]["Recall"]))
