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

def get_truth_dict(file1, file2, window, chroms, min_size=None, max_size=None):
    base_tree = defaultdict(IntervalTree)
    truth_dict = defaultdict(IntervalTree)
    count = 0
    with open(file1) as in_truth:
        for line in in_truth:
            row = line.strip().split('\t')
            #count += 1
            #if count % 10000 == 0:
                #print(count)
            if "mapq<20" in line: 
                continue
            if len(row[7]) < min_size:
                continue
            chrom = row[0]
            if chrom not in chroms:
                continue
            start = int(row[1])
            end = int(row[2])
            key = start
            nearby = base_tree[chrom][start:end]
            nearby_truth = truth_dict[chrom][start:end]
            if len(nearby_truth) == 0:
                count += 1
                truth_dict[chrom][(start-window):(end+window)] = key
    with open(file2) as in_truth:
        for line in in_truth:
            row = line.strip().split('\t')
            #count += 1
            #if count % 10000 == 0:
                #print(count)
            if "mapq<20" in line: 
                continue
            if len(row[7]) < min_size:
                continue
            chrom = row[0]
            if chrom not in chroms:
                continue
            start = int(row[1])
            end = int(row[2])
            key = start
            nearby = base_tree[chrom][start:end]
            nearby_truth = truth_dict[chrom][start:end]
            if len(nearby_truth) == 0:
                count += 1
                truth_dict[chrom][(start-window):(end+window)] = key
    #print(count)
    return truth_dict

def get_vcf_stats(file, truth_dict, min_size):
    seen = {}
    tp = 0
    fp = {}
    fn = 0
    fp_translocations = 0
    for chrom in truth_dict:
        seen[chrom] = defaultdict(int)
    bcf_in = pysam.VariantFile(file)
    for rec in bcf_in.fetch():
        if rec.info["SVTYPE"] == "BND":
            fp_translocations += 1
        elif rec.info["SVTYPE"] == "INS":
            if len(rec.alts[0]) < min_size:
                continue
            chrom = rec.contig
            start = rec.pos
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
    return tp, len(fp), fn, fp_translocations


def get_vcf_distance(file, truth_dict, min_size):
    seen = {}
    for chrom in truth_dict:
        seen[chrom] = defaultdict(int)
    bcf_in = pysam.VariantFile(file)
    for rec in bcf_in.fetch():
        if rec.info["SVTYPE"] == "INS":
            if len(rec.alts[0]) < min_size:
                continue
            chrom = rec.contig
            start = rec.pos
            nearby = truth_dict[chrom][start:start+1]
            closest = 999999999
            pos = -1
            if len(nearby) > 0:
                # Get the closest truth pos for this insert, make sure it only counts towards one
                for item in nearby:
                    if abs(int(item.data)-start) < closest:
                        pos = int(item.data)
                if pos > -1:
                    #print(item.data)
                    if pos in seen[chrom]:
                        seen[chrom][pos] = min(abs(pos-start),seen[chrom][pos])
                    else:
                        seen[chrom][pos] = abs(pos-start)
    return seen

def get_distance_diff(seen_before, seen_after):
    ret_before = []
    ret_after = []
    for chrom in seen_before:
        for pos in seen_before[chrom]:
            if pos in seen_after[chrom]:
                ret_before.append(str(seen_before[chrom][pos]))
                ret_after.append(str(seen_after[chrom][pos]))
    return ",".join(ret_before), ",".join(ret_after)
                


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

parser = argparse.ArgumentParser( description='Get metrics for Sniffles and CuteSV before and after runs')
parser.add_argument('--truth-mat', required=True)
parser.add_argument('--truth-pat', required=True)
parser.add_argument('--sniffles-before', required=True)
parser.add_argument('--sniffles-after', required=True)
parser.add_argument('--cutesv-before', required=True)
parser.add_argument('--cutesv-after', required=True)
parser.add_argument('--window-size', type=int, default=1000)
parser.add_argument('--fraction', type=float, required=True)
parser.add_argument('--rep', type=int, required=True)
parser.add_argument('--chrom-list', default="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY")
parser.add_argument('--insert-stats', required=True)
parser.add_argument('--rt-insert-stats', required=True)
parser.add_argument('--translocation-stats', required=True)
args = parser.parse_args()

# Read in truth data for each repeat type
truth_large = get_truth_dict(args.truth_mat, args.truth_pat, args.window_size, args.chrom_list.split(','), 500)
truth_base = get_truth_dict(args.truth_mat, args.truth_pat, args.window_size, args.chrom_list.split(','), 50)

with open(args.insert_stats, 'w') as out_is, open(args.rt_insert_stats, 'w') as out_rt, open(args.translocation_stats, 'w') as out_tra:
    tp_before, fp_before, fn_before, fp_translocations_before = get_vcf_stats(args.sniffles_before, truth_base, 50)
    tp_after, fp_after, fn_after, fp_translocations_after = get_vcf_stats(args.sniffles_after, truth_base, 50)
    seen_before = get_vcf_distance(args.sniffles_before, truth_base, 50)
    seen_after = get_vcf_distance(args.sniffles_after, truth_base, 50)
    diff_before_sniffles_50, diff_after_sniffles_50 = get_distance_diff(seen_before, seen_after)
    sniffles_before_pb_precision = tp_before/(tp_before+fp_before)
    sniffles_after_pb_precision = tp_after/(tp_after+fp_after)
    sniffles_before_pb_recall = tp_before/(tp_before+fn_before)
    sniffles_after_pb_recall = tp_after/(tp_after+fn_after)

    out_tra.write(str(args.fraction)+"\t"+str(args.rep)+"\t"+str(fp_translocations_before)+"\t"+str(fp_translocations_after)+"\t")
    #print(str(args.fraction)+"\t"+str(args.rep)+"\t"+str(fp_translocations_before)+"\t"+str(fp_translocations_after))
    out_is.write(str(args.fraction)+"\t"+str(args.rep)+"\t"+str(tp_before)+"\t"+str(fp_before)+"\t"+str(fn_before)+"\t"+str(sniffles_before_pb_precision)+"\t"+str(sniffles_before_pb_recall)+"\t"+str(tp_after)+"\t"+str(fp_after)+"\t"+str(fn_after)+"\t"+str(sniffles_after_pb_precision)+"\t"+str(sniffles_after_pb_recall)+"\t")
    #print(str(args.fraction)+"\t"+str(args.rep)+"\t"+str(tp_before)+"\t"+str(fp_before)+"\t"+str(fn_before)+"\t"+str(sniffles_before_pb_precision)+"\t"+str(sniffles_before_pb_recall)+"\t"+str(tp_after)+"\t"+str(fp_after)+"\t"+str(fn_after)+"\t"+str(sniffles_after_pb_precision)+"\t"+str(sniffles_after_pb_recall))

    tp_before, fp_before, fn_before, fp_translocations_before = get_vcf_stats(args.sniffles_before, truth_large, 500)
    tp_after, fp_after, fn_after, fp_translocations_after = get_vcf_stats(args.sniffles_after, truth_large, 500)
    seen_before = get_vcf_distance(args.sniffles_before, truth_large, 500)
    seen_after = get_vcf_distance(args.sniffles_after, truth_large, 500)
    diff_before_sniffles_500, diff_after_sniffles_500 = get_distance_diff(seen_before, seen_after)
    sniffles_before_pb_precision = tp_before/(tp_before+fp_before)
    sniffles_after_pb_precision = tp_after/(tp_after+fp_after)
    sniffles_before_pb_recall = tp_before/(tp_before+fn_before)
    sniffles_after_pb_recall = tp_after/(tp_after+fn_after)

    out_rt.write(str(args.fraction)+"\t"+str(args.rep)+"\t"+str(tp_before)+"\t"+str(fp_before)+"\t"+str(fn_before)+"\t"+str(sniffles_before_pb_precision)+"\t"+str(sniffles_before_pb_recall)+"\t"+str(tp_after)+"\t"+str(fp_after)+"\t"+str(fn_after)+"\t"+str(sniffles_after_pb_precision)+"\t"+str(sniffles_after_pb_recall)+"\t")
    #print(str(args.fraction)+"\t"+str(args.rep)+"\t"+str(tp_before)+"\t"+str(fp_before)+"\t"+str(fn_before)+"\t"+str(sniffles_before_pb_precision)+"\t"+str(sniffles_before_pb_recall)+"\t"+str(tp_after)+"\t"+str(fp_after)+"\t"+str(fn_after)+"\t"+str(sniffles_after_pb_precision)+"\t"+str(sniffles_after_pb_recall))

    tp_before, fp_before, fn_before, fp_translocations_before = get_vcf_stats(args.cutesv_before, truth_base, 50)
    tp_after, fp_after, fn_after, fp_translocations_after = get_vcf_stats(args.cutesv_after, truth_base, 50)
    seen_before = get_vcf_distance(args.cutesv_before, truth_base, 50)
    seen_after = get_vcf_distance(args.cutesv_after, truth_base, 50)
    diff_before_cutesv_50, diff_after_cutesv_50 = get_distance_diff(seen_before, seen_after)
    cutesv_before_pb_precision = tp_before/(tp_before+fp_before)
    cutesv_after_pb_precision = tp_after/(tp_after+fp_after)
    cutesv_before_pb_recall = tp_before/(tp_before+fn_before)
    cutesv_after_pb_recall = tp_after/(tp_after+fn_after)

    out_is.write(str(tp_before)+"\t"+str(fp_before)+"\t"+str(fn_before)+"\t"+str(cutesv_before_pb_precision)+"\t"+str(cutesv_before_pb_recall)+"\t"+str(tp_after)+"\t"+str(fp_after)+"\t"+str(fn_after)+"\t"+str(cutesv_after_pb_precision)+"\t"+str(cutesv_after_pb_recall)+"\t"+diff_before_sniffles_50+"\t"+diff_after_sniffles_50+"\t"+diff_before_cutesv_50+"\t"+diff_after_cutesv_50+"\n")
    #print(str(args.fraction)+"\t"+str(args.rep)+"\t"+str(tp_before)+"\t"+str(fp_before)+"\t"+str(fn_before)+"\t"+str(sniffles_before_pb_precision)+"\t"+str(sniffles_before_pb_recall)+"\t"+str(tp_after)+"\t"+str(fp_after)+"\t"+str(fn_after)+"\t"+str(sniffles_after_pb_precision)+"\t"+str(sniffles_after_pb_recall))

    tp_before, fp_before, fn_before, fp_translocations_before = get_vcf_stats(args.cutesv_before, truth_large, 500)
    tp_after, fp_after, fn_after, fp_translocations_after = get_vcf_stats(args.cutesv_after, truth_large, 500)
    seen_before = get_vcf_distance(args.cutesv_before, truth_large, 500)
    seen_after = get_vcf_distance(args.cutesv_after, truth_large, 500)
    diff_before_cutesv_500, diff_after_cutesv_500 = get_distance_diff(seen_before, seen_after)
    cutesv_before_pb_precision = tp_before/(tp_before+fp_before)
    cutesv_after_pb_precision = tp_after/(tp_after+fp_after)
    cutesv_before_pb_recall = tp_before/(tp_before+fn_before)
    cutesv_after_pb_recall = tp_after/(tp_after+fn_after)

    out_tra.write(str(fp_translocations_before)+"\t"+str(fp_translocations_after)+"\n")
    #print(str(fp_translocations_before)+"\t"+str(fp_translocations_after))
    out_rt.write(str(tp_before)+"\t"+str(fp_before)+"\t"+str(fn_before)+"\t"+str(cutesv_before_pb_precision)+"\t"+str(cutesv_before_pb_recall)+"\t"+str(tp_after)+"\t"+str(fp_after)+"\t"+str(fn_after)+"\t"+str(cutesv_after_pb_precision)+"\t"+str(cutesv_after_pb_recall)+"\t"+diff_before_sniffles_500+"\t"+diff_after_sniffles_500+"\t"+diff_before_cutesv_500+"\t"+diff_after_cutesv_500+"\n")
    #print(str(args.fraction)+"\t"+str(args.rep)+"\t"+str(tp_before)+"\t"+str(fp_before)+"\t"+str(fn_before)+"\t"+str(sniffles_before_pb_precision)+"\t"+str(sniffles_before_pb_recall)+"\t"+str(tp_after)+"\t"+str(fp_after)+"\t"+str(fn_after)+"\t"+str(sniffles_after_pb_precision)+"\t"+str(sniffles_after_pb_recall))

