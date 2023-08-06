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
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

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

def get_truth_dict(file, window, chroms, min_size=None, max_size=None):
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

def plot_tra_bar(name, x, y1, y2, caller):
    X_axis = np.arange(len(x))
    plt.bar(X_axis - 0.2, y1, 0.4, label = 'Before')
    plt.bar(X_axis + 0.2, y2, 0.4, label = 'After')
    plt.xticks(X_axis, x)
    plt.xlabel("Coverage")
    plt.ylabel("False Positive Translocations")
    plt.title(caller+" FP Translocations before and after ReAlignment")
    plt.legend()
    fig = plt.gcf()
    fig.savefig(name)
    plt.clf()

def plot_precision_recall_bar(name, x, y1, y2, ylabel, title):
    X_axis = np.arange(len(x))
    plt.bar(X_axis - 0.2, y1, 0.4, label = 'Before')
    plt.bar(X_axis + 0.2, y2, 0.4, label = 'After')
    plt.xticks(X_axis, x)
    plt.xlabel("Coverage")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    fig = plt.gcf()
    fig.savefig(name)
    plt.clf()

def plot_pr_gain(name, x, y1, ylabel, title):
    X_axis = np.arange(len(x))
    plt.bar(X_axis - 0.2, y1*100, 0.4, label = 'Gain')
    plt.xticks(X_axis, x)
    plt.xlabel("Coverage")
    plt.ylabel(ylabel)
    plt.title(title)
    fig = plt.gcf()
    fig.savefig(name)
    plt.clf()

def get_tp_fp_fn(prefix, coverage, rep, realigned, vcf, tp_svs, chroms, min_size=None, max_size=None):
    tp = 0
    fp = 0
    fn = 0
    seen = {}
    add = ""
    if realigned:
        add = "Realigned_"
    if not vcf:
        # BEDPE
        my_file = Path(prefix+"_"+add+"HG002_"+coverage+"_"+rep+"/"+add+"HG002_"+coverage+"_"+rep+".bedpe")
        if my_file.is_file():
            with open(prefix+"_"+add+"HG002_"+coverage+"_"+rep+"/"+add+"HG002_"+coverage+"_"+rep+".bedpe", 'r') as in_bedpe:
                count = 0
                for line in in_bedpe:
                    if count == 0:
                        count += 1
                        continue
                    row = line.strip().split('\t')
                    if len(row) < 16:
                        continue
                    if row[10] == "INS":
                        if row[16] == "NA":
                            continue
                        if min_size is not None:
                            if int(row[16]) < min_size:
                                continue
                        if max_size is not None:
                            if int(row[16]) > max_size:
                                continue
                        nearby = tp_svs[row[0]][int(row[1]):int(row[1])+1]
                        if len(nearby) > 0:
                            #TP
                            for item in nearby:
                                if row[0] not in seen:
                                    seen[row[0]] = {}
                                seen[row[0]][item.data] = 1
                        else:
                            #FP
                            fp += 1
        for chrom in chroms:
            if chrom not in seen:
                for pos in tp_svs[chrom]:
                    fn += 1
            else:
                for pos in tp_svs[chrom]:
                    if pos.data in seen[chrom]:
                        tp += 1
                    else:
                        fn += 1
    else:
        my_file = Path(prefix+"_"+add+"HG002_"+coverage+"_"+rep+"/"+add+"HG002_"+coverage+"_"+rep+".vcf")
        if my_file.is_file():
            vcf_in = pysam.VariantFile(prefix+"_"+add+"HG002_"+coverage+"_"+rep+"/"+add+"HG002_"+coverage+"_"+rep+".vcf")
            for rec in vcf_in.fetch():
                if rec.info["SVTYPE"] == "INS":
                    if min_size is not None:
                        if len(rec.alts[0]) < min_size:
                            continue
                    if max_size is not None:
                        if len(rec.alts[0]) > max_size:
                            continue
                    nearby = tp_svs[rec.contig][rec.pos:rec.pos+1]
                    if len(nearby) > 0:
                        #TP
                        for item in nearby:
                            if rec.contig not in seen:
                                seen[rec.contig] = {}
                            seen[rec.contig][item.data] = 1
                    else:
                        #FP
                        fp += 1
        for chrom in chroms:
            if chrom not in seen:
                for pos in tp_svs[chrom]:
                    fn += 1
            else:
                for pos in tp_svs[chrom]:
                    if pos.data in seen[chrom]:
                        tp += 1
                    else:
                        fn += 1
    return tp, fp, fn

def get_tra(prefix, coverage, rep, realigned, vcf):
    fp_count = 0
    add = ""
    if realigned:
        add = "Realigned_"
    if not vcf:
        # BEDPE
        my_file = Path(prefix+"_"+add+"HG002_"+coverage+"_"+rep+"/"+add+"HG002_"+coverage+"_"+rep+".bedpe")
        if my_file.is_file():
            with open(prefix+"_"+add+"HG002_"+coverage+"_"+rep+"/"+add+"HG002_"+coverage+"_"+rep+".bedpe", 'r') as in_bedpe:
                count = 0
                for line in in_bedpe:
                    if count == 0:
                        count += 1
                        continue
                    row = line.strip().split('\t')
                    if len(row) < 10:
                        continue
                    if row[10] == "TRA":
                        fp_count += 1
    else:
        my_file = Path(prefix+"_"+add+"HG002_"+coverage+"_"+rep+"/"+add+"HG002_"+coverage+"_"+rep+".vcf")
        if my_file.is_file():
            vcf_in = pysam.VariantFile(prefix+"_"+add+"HG002_"+coverage+"_"+rep+"/"+add+"HG002_"+coverage+"_"+rep+".vcf")
            for rec in vcf_in.fetch():
                if rec.info["SVTYPE"] == "BND":
                    fp_count += 1
    return fp_count

parser = argparse.ArgumentParser( description='Get metrics for xtea and tldr before and after runs')
#parser.add_argument('--truth-line', required=True)
parser.add_argument('--sniffles', default="Sniffles")
parser.add_argument('--cutesv', default="CuteSV")
parser.add_argument('--truth', required=True)
parser.add_argument('--window-size', type=int, default=1000)
parser.add_argument('--chrom-list', default="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY")
args = parser.parse_args()

# First get FP Translocation calls before and after
sniffles_before = []
sniffles_after = []
cute_sv_before = []
cute_sv_after = []
for coverage in ["0.05", "0.1", "0.15", "0.2", "0.25", "0.3", "0.35", "0.4", "0.45", "0.5", "0.55", "0.6", "0.65", "0.7", "0.75", "0.8", "0.85", "0.9", "0.95"]:
    sniffles_before_tra_fp = []
    sniffles_after_tra_fp = []
    cute_sv_before_tra_fp = []
    cute_sv_after_tra_fp = []
    for rep in ["1", "2", "3"]:
        sniffles_before_tra_fp.append(get_tra(args.sniffles, coverage, rep, False, True))
        sniffles_after_tra_fp.append(get_tra(args.sniffles, coverage, rep, True, True))
        cute_sv_before_tra_fp.append(get_tra(args.cutesv, coverage, rep, False, True))
        cute_sv_after_tra_fp.append(get_tra(args.cutesv, coverage, rep, True, True))
    sniffles_before.append(np.mean(sniffles_before_tra_fp))
    sniffles_after.append(np.mean(sniffles_after_tra_fp))
    cute_sv_before.append(np.mean(cute_sv_before_tra_fp))
    cute_sv_after.append(np.mean(cute_sv_after_tra_fp))


# Plot Before and after for each repeat family
print("Translocations")
x = np.array(["3", "6", "9", "12", "15", "18", "21", "24", "27", "30", "33", "36", "39", "42", "45", "48", "51", "54", "57"])
plot_tra_bar('/u/adcosta/transfer/sniffles_tra.png', x, np.array(sniffles_before), np.array(sniffles_after), "Sniffles")
plot_tra_bar('/u/adcosta/transfer/cutesv_tra.png', x, np.array(cute_sv_before), np.array(cute_sv_after), "CuteSV")


chroms = {}
for chrom in args.chrom_list.split(','):
    chroms[chrom] = 1

for pair in [[None, None], [100, 500], [500, 2500], [2500, None]]:
    print(pair)
    min_size = 30
    max_size = None
    suffix = ""
    if pair[0] is not None and pair[1] is not None:
        suffix = "_"+str(pair[0])+"_"+str(pair[1])
    elif pair[0] is not None:
        suffix = "_"+str(pair[0])+"_Plus"
    if pair[0] is not None:
        min_size = pair[0]
    if pair[1] is not None:
        max_size = pair[1]
    # Read in a list of SV positions for HG002
    tp_svs = get_truth_dict(args.truth, args.window_size, chroms, min_size, max_size)
    # Use list as TP set for Sniffles and CuteSV calls, Compare Preceision and Recall before and after
    sniffles_before_precision = []
    sniffles_after_precision = []
    cute_sv_before_precision = []
    cute_sv_after_precision = []
    sniffles_before_recall = []
    sniffles_after_recall = []
    cute_sv_before_recall = []
    cute_sv_after_recall = []
    cute_sv_precision_gain = []
    sniffles_precision_gain = []
    cute_sv_recall_gain = []
    sniffles_recall_gain = []
    for coverage in ["0.05", "0.1", "0.15", "0.2", "0.25", "0.3", "0.35", "0.4", "0.45", "0.5", "0.55", "0.6", "0.65", "0.7", "0.75", "0.8", "0.85", "0.9", "0.95"]:
        #print(coverage)
        sniffles_before_cov_precision = []
        sniffles_after_cov_precision = []
        cute_sv_before_cov_precision = []
        cute_sv_after_cov_precision = []
        sniffles_before_cov_recall = []
        sniffles_after_cov_recall = []
        cute_sv_before_cov_recall = []
        cute_sv_after_cov_recall = []
        for rep in ["1", "2", "3"]:
            #print(rep)
            tp, fp, fn = get_tp_fp_fn(args.sniffles, coverage, rep, False, True, tp_svs, chroms, min_size, max_size)
            if tp == 0 and fp == 0:
                sniffles_before_cov_precision.append(0)
            else:
                sniffles_before_cov_precision.append(tp/(tp+fp))
            if tp == 0 and fn == 0:
                sniffles_before_cov_recall.append(0)
            else:
                sniffles_before_cov_recall.append(tp/(tp+fn))
            tp, fp, fn = get_tp_fp_fn(args.sniffles, coverage, rep, True, True, tp_svs, chroms, min_size, max_size)
            if tp == 0 and fp == 0:
                sniffles_after_cov_precision.append(0)
            else:
                sniffles_after_cov_precision.append(tp/(tp+fp))
            if tp == 0 and fn == 0:
                sniffles_after_cov_recall.append(0)
            else:
                sniffles_after_cov_recall.append(tp/(tp+fn))
            tp, fp, fn = get_tp_fp_fn(args.cutesv, coverage, rep, False, True, tp_svs, chroms, min_size, max_size)
            if tp == 0 and fp == 0:
                cute_sv_before_cov_precision.append(0)
            else:
                cute_sv_before_cov_precision.append(tp/(tp+fp))
            if tp == 0 and fn == 0:
                cute_sv_before_cov_recall.append(0)
            else:
                cute_sv_before_cov_recall.append(tp/(tp+fn))
            tp, fp, fn = get_tp_fp_fn(args.cutesv, coverage, rep, True, True, tp_svs, chroms, min_size, max_size)
            if tp == 0 and fp == 0:
                cute_sv_after_cov_precision.append(0)
            else:
                cute_sv_after_cov_precision.append(tp/(tp+fp))
            if tp == 0 and fn == 0:
                cute_sv_after_cov_recall.append(0)
            else:
                cute_sv_after_cov_recall.append(tp/(tp+fn))
        sniffles_before_precision.append(np.mean(sniffles_before_cov_precision))
        sniffles_after_precision.append(np.mean(sniffles_after_cov_precision))
        cute_sv_before_precision.append(np.mean(cute_sv_before_cov_precision))
        cute_sv_after_precision.append(np.mean(cute_sv_after_cov_precision))
        sniffles_before_recall.append(np.mean(sniffles_before_cov_recall))
        sniffles_after_recall.append(np.mean(sniffles_after_cov_recall))
        cute_sv_before_recall.append(np.mean(cute_sv_before_cov_recall))
        cute_sv_after_recall.append(np.mean(cute_sv_after_cov_recall))
        cute_sv_precision_gain.append(np.mean(np.array(cute_sv_after_cov_precision) - np.array(cute_sv_before_cov_precision)))
        sniffles_precision_gain.append(np.mean(np.array(sniffles_after_cov_precision) - np.array(sniffles_before_cov_precision)))
        cute_sv_recall_gain.append(np.mean(np.array(cute_sv_after_cov_recall) - np.array(cute_sv_before_cov_recall)))
        sniffles_recall_gain.append(np.mean(np.array(sniffles_after_cov_recall) - np.array(sniffles_before_cov_recall)))
    print("Precision")
    plot_precision_recall_bar('/u/adcosta/transfer/plots/sniffles_precision'+suffix+'.png', x, np.array(sniffles_before_precision), np.array(sniffles_after_precision), "Precision", "Sniffles insert precision before and after ReAlignment")
    plot_precision_recall_bar('/u/adcosta/transfer/plots/cutesv_precision'+suffix+'.png', x, np.array(cute_sv_before_precision), np.array(cute_sv_after_precision), "Precision", "CuteSV insert precision before and after ReAlignment")
    print("Recall")
    plot_precision_recall_bar('/u/adcosta/transfer/plots/sniffles_recall'+suffix+'.png', x, np.array(sniffles_before_recall), np.array(sniffles_after_recall), "Recall", "Sniffles insert recall before and after ReAlignment")
    plot_precision_recall_bar('/u/adcosta/transfer/plots/cutesv_recall'+suffix+'.png', x, np.array(cute_sv_before_recall), np.array(cute_sv_after_recall), "Recall", "CuteSV insert recall before and after ReAlignment")
    print("Recall Gain")
    plot_pr_gain('/u/adcosta/transfer/plots/sniffles_recall_gain'+suffix+'.png', x, np.array(sniffles_recall_gain), "Recall % Gain", "Sniffles insert recall gain after ReAlignment")
    plot_pr_gain('/u/adcosta/transfer/plots/cutesv_recall_gain'+suffix+'.png', x, np.array(cute_sv_recall_gain), "Recall % Gain", "CuteSV insert recall gain after ReAlignment")
    print("Precision Gain")
    plot_pr_gain('/u/adcosta/transfer/plots/sniffles_precision_gain'+suffix+'.png', x, np.array(sniffles_precision_gain), "Precision % Gain", "Sniffles insert precision gain after ReAlignment")
    plot_pr_gain('/u/adcosta/transfer/plots/cutesv_precision_gain'+suffix+'.png', x, np.array(cute_sv_precision_gain), "Precision % Gain", "CuteSV insert precision gain after ReAlignment")

