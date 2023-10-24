import pysam
import argparse
import sys
import csv
import os
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

def get_nearby(chrom, start, end, seen):
    if chrom in seen:
        count = 0
        i = start
        while i < end:
            if i in seen[chrom]:
                count += seen[chrom][i]
            i += 1
        if count >= 3:
            return True
        return False
    return True

def get_nearby_tp(chrom, start, end, seen):
    if chrom in seen:
        count = 0
        i = start
        while i < end:
            if i in seen[chrom]:
                count += seen[chrom][i]
            i += 1
        if count > 0:
            return True
        return False
    return True

def get_truth_dict(file1, window, chroms, min_size=None, max_size=None):
    truth_all = defaultdict(IntervalTree)
    truth_rt = defaultdict(IntervalTree)
    count = 0
    with open(file1) as in_truth:
        for line in in_truth:
            if count == 0:
                count += 1
                continue
            if "PASS" not in line: 
                continue
            row = line.strip().split('\t')
            chrom = row[0]
            if chrom not in chroms:
                continue
            start = int(row[1])-window
            end = int(row[2])+window
            key = start
            annotation = "NoAnnotation"
            if len(row) > 9:
                if len(row[9]) > 0:
                    annotation = row[9]
            truth_all[chrom][start:end] = [int(row[1]), int(row[2]), annotation]
            if "alu" in annotation or "line" in annotation or "sva" in annotation:
                truth_rt[chrom][start:end] = [int(row[1]), int(row[2]), annotation]
    return truth_all, truth_rt

def get_family(item1, item2):
    if "alu" in item1.lower() and "alu" in item2.lower():
        return True
    if "sva" in item1.lower() and "sva" in item2.lower():
        return True
    if "line" in item1.lower() and "line" in item2.lower():
        return True
    return False

def get_xtea_stats_rt(file_prefix, truth_dict):
    seen = {}
    tp = 0
    fp = {}
    fn = 0
    tp_distance = 0
    for chrom in truth_dict:
        seen[chrom] = defaultdict(int)
    if not os.path.isfile(file_prefix+".merged_ALU.txt"):
        return tp, fp, fn
    with open(file_prefix+".merged_ALU.txt") as in_calls:
        for line in in_calls:
            row = line.strip().split('\t')
            chrom = row[0]
            start = int(row[1])
            # Check if the position is within window size of a truth call
            nearby = truth_dict[chrom][start:start+1]
            if len(nearby) > 0:
                # Have a tp
                found = False
                for item in nearby:
                    #print(item.data)
                    if get_family(item.data[2], row[2]):
                        found = True
                        if int(item.data[0]) in seen[chrom]:
                            seen[chrom][int(item.data[0])] = min(abs(int(item.data[0])-start),seen[chrom][int(item.data[0])])
                        else:
                            seen[chrom][int(item.data[0])] = abs(int(item.data[0])-start)
                if not found:
                    # Insertion was nearby to a TP call but did not have the same repeat family 
                    if chrom in truth_dict:
                        fp[chrom+":"+str(start)] = 1
                        print(chrom+"\t"+str(start)+"\t"+str(start+1)+"\tFP\txTea-Long")
            else:
                # Have a fp
                if chrom in truth_dict:
                    fp[chrom+":"+str(start)] = 1
                    print(chrom+"\t"+str(start)+"\t"+str(start+1)+"\tFP\txTea-Long")
    with open(file_prefix+".merged_LINE1.txt") as in_calls:
        for line in in_calls:
            row = line.strip().split('\t')
            chrom = row[0]
            start = int(row[1])
            # Check if the position is within window size of a truth call
            nearby = truth_dict[chrom][start:start+1]
            if len(nearby) > 0:
                # Have a tp
                found = False
                for item in nearby:
                    #print(item.data)
                    if get_family(item.data[2], row[2]):
                        found = True
                        if int(item.data[0]) in seen[chrom]:
                            seen[chrom][int(item.data[0])] = min(abs(int(item.data[0])-start),seen[chrom][int(item.data[0])])
                        else:
                            seen[chrom][int(item.data[0])] = abs(int(item.data[0])-start)
                if not found:
                    # Insertion was nearby to a TP call but did not have the same repeat family 
                    if chrom in truth_dict:
                        fp[chrom+":"+str(start)] = 1
                        print(chrom+"\t"+str(start)+"\t"+str(start+1)+"\tFP\txTea-Long")
            else:
                # Have a fp
                if chrom in truth_dict:
                    fp[chrom+":"+str(start)] = 1
                    print(chrom+"\t"+str(start)+"\t"+str(start+1)+"\tFP\txTea-Long")
    with open(file_prefix+".merged_HERV.txt") as in_calls:
        for line in in_calls:
            row = line.strip().split('\t')
            chrom = row[0]
            start = int(row[1])
            # Check if the position is within window size of a truth call
            nearby = truth_dict[chrom][start:start+1]
            if len(nearby) > 0:
                # Have a tp
                found = False
                for item in nearby:
                    #print(item.data)
                    if get_family(item.data[2], row[2]):
                        found = True
                        if int(item.data[0]) in seen[chrom]:
                            seen[chrom][int(item.data[0])] = min(abs(int(item.data[0])-start),seen[chrom][int(item.data[0])])
                        else:
                            seen[chrom][int(item.data[0])] = abs(int(item.data[0])-start)
                if not found:
                    # Insertion was nearby to a TP call but did not have the same repeat family 
                    if chrom in truth_dict:
                        fp[chrom+":"+str(start)] = 1
                        print(chrom+"\t"+str(start)+"\t"+str(start+1)+"\tFP\txTea-Long")
            else:
                # Have a fp
                if chrom in truth_dict:
                    fp[chrom+":"+str(start)] = 1
                    print(chrom+"\t"+str(start)+"\t"+str(start+1)+"\tFP\txTea-Long")
    with open(file_prefix+".merged_SVA.txt") as in_calls:
        for line in in_calls:
            row = line.strip().split('\t')
            chrom = row[0]
            start = int(row[1])
            # Check if the position is within window size of a truth call
            nearby = truth_dict[chrom][start:start+1]
            if len(nearby) > 0:
                # Have a tp
                found = False
                for item in nearby:
                    #print(item.data)
                    if get_family(item.data[2], row[2]):
                        found = True
                        if int(item.data[0]) in seen[chrom]:
                            seen[chrom][int(item.data[0])] = min(abs(int(item.data[0])-start),seen[chrom][int(item.data[0])])
                        else:
                            seen[chrom][int(item.data[0])] = abs(int(item.data[0])-start)
                if not found:
                    # Insertion was nearby to a TP call but did not have the same repeat family 
                    if chrom in truth_dict:
                        fp[chrom+":"+str(start)] = 1
                        print(chrom+"\t"+str(start)+"\t"+str(start+1)+"\tFP\txTea-Long")
            else:
                # Have a fp
                if chrom in truth_dict:
                    fp[chrom+":"+str(start)] = 1  
                    print(chrom+"\t"+str(start)+"\t"+str(start+1)+"\tFP\txTea-Long")  
    for chrom in truth_dict:
        for item in truth_dict[chrom]:
            if int(item.data[0]) not in seen[chrom]:
                print(chrom+"\t"+str(item.data[0])+"\t"+str(item.data[0]+1)+"\tFN\txTea-Long")
                fn += 1
    for chrom in seen:
        for pos in seen[chrom]:
            tp += 1
            print(chrom+"\t"+str(pos)+"\t"+str(pos+1)+"\tTP\txTea-Long")
            tp_distance += seen[chrom][pos]
    return tp, fp, fn


def get_tldr_stats_rt(file, truth_dict):
    seen = {}
    tp = 0
    fp = 1
    fn = 0
    for chrom in truth_dict:
        seen[chrom] = defaultdict(int)
    count = 0
    if not os.path.isfile(file):
        return tp, fp, fn
    with open(file) as in_calls:
        for line in in_calls:
            if count == 0:
                count += 1
                continue
            if "PASS" not in line:
                # tldr has not annotated this insert to a repeat family or annotated to less than 50%
                continue
            row = line.strip().split('\t')
            chrom = row[1]
            start = int(row[2])
            end = int(row[3])
            # Check if the position is within window size of a truth call
            nearby = truth_dict[chrom][start:end]
            if len(nearby) > 0:
                # Have a tp
                found = False
                for item in nearby:
                    #print(item.data)
                    if get_family(item.data[2], row[5]):
                        found = True
                        seen[chrom][int(item.data[0])] = 1
                if not found:
                    # Insertion was nearby to a TP call but did not have the same repeat family 
                    if chrom in truth_dict:
                        print(chrom+"\t"+str(start)+"\t"+str(end)+"\tFP\ttldr")
                        fp += 1
            else:
                # Have a fp
                if chrom in truth_dict:
                    print(chrom+"\t"+str(start)+"\t"+str(end)+"\tFP\ttldr")
                    fp += 1
        for chrom in truth_dict:
            for item in truth_dict[chrom]:
                if int(item.data[0]) not in seen[chrom]:
                    print(chrom+"\t"+str(item.data[0])+"\t"+str(item.data[0]+1)+"\tFN\ttldr")
                    fn += 1
        for chrom in seen:
            for pos in seen[chrom]:
                print(chrom+"\t"+str(pos)+"\t"+str(pos+1)+"\tTP\ttldr")
                tp += 1
    return tp, fp, fn


def get_somrit_stats_rt(file, window, truth_dict):
    base_tree = defaultdict(IntervalTree)
    fp_dict = defaultdict(IntervalTree)
    seen = {}
    all_inserts = {}
    tp = 0
    fp = 0
    fn = 0
    for chrom in truth_dict:
        seen[chrom] = defaultdict(int)
        all_inserts[chrom] = defaultdict(int)
    count = 0
    with open(file) as in_calls:
        for line in in_calls:
            if count == 0:
                count += 1
                continue
            if "PASS" not in line:
                continue
            row = line.strip().split('\t')
            chrom = row[0]
            start = int(row[1])
            end = int(row[2])
            # Check if the position is within window size of a truth call
            nearby = truth_dict[chrom][start:end]
            if len(nearby) > 0:
                # Have a tp
                for item in nearby:
                    #print(item.data)
                    if get_family(item.data[2], row[7]):
                        seen[chrom][int(item.data[0])] = 1
            if chrom in truth_dict:
                # Have a fp
                all_inserts[chrom][start] += len(row[5].split(',')) 
    count = 0
    with open(file) as in_calls:
        for line in in_calls:
            if count == 0:
                count += 1
                continue
            if "PASS" not in line:
                continue
            row = line.strip().split('\t')
            chrom = row[0]
            start = int(row[1])
            end = int(row[2])
            nearby = truth_dict[chrom][start:end]
            found = False
            if len(nearby) > 0:
                # Have a tp
                for item in nearby:
                    #print(item.data)
                    if get_family(item.data[2], row[7]):
                        found = True
            if found:
                continue
            nearby = get_nearby(chrom, start-50, end+50, all_inserts)
            nearby_fp = fp_dict[chrom][start:end]
            if nearby and len(nearby_fp) == 0:
                fp_dict[chrom][(start-50):(end+50)] = start
    for chrom in truth_dict:
        for item in truth_dict[chrom]:
            if int(item.data[0]) not in seen[chrom]:
                print(chrom+"\t"+str(item.data[0])+"\t"+str(item.data[0]+1)+"\tFN\tsomrit")
                fn += 1
    for chrom in seen:
        for pos in seen[chrom]:
            print(chrom+"\t"+str(pos)+"\t"+str(pos+1)+"\tTP\tsomrit")
            tp += 1
    for chrom in truth_dict:
        for pos in fp_dict[chrom]:
            print(chrom+"\t"+str(pos.data)+"\t"+str(pos.data+1)+"\tFP\tsomrit")
            fp += 1
    return tp, fp, fn



def get_xtea_stats_all(file_prefix, truth_dict):
    seen = {}
    tp = 0
    fp = {}
    fn = 0
    tp_distance = 0
    for chrom in truth_dict:
        seen[chrom] = defaultdict(int)
    if not os.path.isfile(file_prefix+".merged_ALU.txt"):
        return tp, fp, fn
    with open(file_prefix+".merged_ALU.txt") as in_calls:
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
                    if int(item.data[0]) in seen[chrom]:
                        seen[chrom][int(item.data[0])] = min(abs(int(item.data[0])-start),seen[chrom][int(item.data[0])])
                    else:
                        seen[chrom][int(item.data[0])] = abs(int(item.data[0])-start)
            else:
                # Have a fp
                if chrom in truth_dict:
                    fp[chrom+":"+str(start)] = 1
    with open(file_prefix+".merged_LINE1.txt") as in_calls:
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
                    if int(item.data[0]) in seen[chrom]:
                        seen[chrom][int(item.data[0])] = min(abs(int(item.data[0])-start),seen[chrom][int(item.data[0])])
                    else:
                        seen[chrom][int(item.data[0])] = abs(int(item.data[0])-start)
            else:
                # Have a fp
                if chrom in truth_dict:
                    fp[chrom+":"+str(start)] = 1
    with open(file_prefix+".merged_HERV.txt") as in_calls:
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
                    if int(item.data[0]) in seen[chrom]:
                        seen[chrom][int(item.data[0])] = min(abs(int(item.data[0])-start),seen[chrom][int(item.data[0])])
                    else:
                        seen[chrom][int(item.data[0])] = abs(int(item.data[0])-start)
            else:
                # Have a fp
                if chrom in truth_dict:
                    fp[chrom+":"+str(start)] = 1
    with open(file_prefix+".merged_SVA.txt") as in_calls:
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
                    if int(item.data[0]) in seen[chrom]:
                        seen[chrom][int(item.data[0])] = min(abs(int(item.data[0])-start),seen[chrom][int(item.data[0])])
                    else:
                        seen[chrom][int(item.data[0])] = abs(int(item.data[0])-start)
            else:
                # Have a fp
                if chrom in truth_dict:
                    fp[chrom+":"+str(start)] = 1
    for chrom in truth_dict:
        for item in truth_dict[chrom]:
            if int(item.data[0]) not in seen[chrom]:
                fn += 1
    for chrom in seen:
        for pos in seen[chrom]:
            tp += 1
            tp_distance += seen[chrom][pos]
    return tp, fp, fn


def get_tldr_stats_all(file, truth_dict):
    seen = {}
    tp = 0
    fp = 1
    fn = 0
    for chrom in truth_dict:
        seen[chrom] = defaultdict(int)
    count = 0
    if not os.path.isfile(file):
        return tp, fp, fn
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
            nearby = truth_dict[chrom][start:end]
            if len(nearby) > 0:
                # Have a tp
                for item in nearby:
                    #print(item.data)
                    seen[chrom][int(item.data[0])] = 1
            else:
                # Have a fp
                if chrom in truth_dict:
                    fp += 1
        for chrom in truth_dict:
            for item in truth_dict[chrom]:
                if int(item.data[0]) not in seen[chrom]:
                    fn += 1
        for chrom in seen:
            for pos in seen[chrom]:
                tp += 1
    return tp, fp, fn


def get_somrit_stats_all(file, window, truth_dict):
    base_tree = defaultdict(IntervalTree)
    fp_dict = defaultdict(IntervalTree)
    seen = {}
    all_inserts = {}
    tp = 0
    fp = 0
    fn = 0
    for chrom in truth_dict:
        seen[chrom] = defaultdict(int)
        all_inserts[chrom] = defaultdict(int)
    count = 0
    with open(file) as in_calls:
        for line in in_calls:
            if count == 0:
                count += 1
                continue
            row = line.strip().split('\t')
            #if "PASS" not in line:
            #    if row[7] != "No_Mapping":
            #        continue
            chrom = row[0]
            start = int(row[1])
            end = int(row[2])
            # Check if the position is within window size of a truth call
            nearby = truth_dict[chrom][start:end]
            if len(nearby) > 0:
                # Have a tp
                for item in nearby:
                    #print(item.data)
                    seen[chrom][int(item.data[0])] = 1
            if chrom in truth_dict:
                all_inserts[chrom][start] += len(row[5].split(',')) 
    count = 0
    with open(file) as in_calls:
        for line in in_calls:
            if count == 0:
                count += 1
                continue
            row = line.strip().split('\t')
            #if "PASS" not in line:
            #    if row[7] != "No_Mapping":
            #        continue
            chrom = row[0]
            start = int(row[1])
            end = int(row[2])
            nearby = get_nearby(chrom, start-50, end+50, all_inserts)
            nearby_truth = truth_dict[chrom][start:end]
            nearby_fp = fp_dict[chrom][start:end]
            if len(nearby_truth) == 0 and nearby and len(nearby_fp) == 0:
                fp_dict[chrom][(start-50):(end+50)] = start
    for chrom in truth_dict:
        for item in truth_dict[chrom]:
            if int(item.data[0]) not in seen[chrom]:
                fn += 1
    for chrom in seen:
        for pos in seen[chrom]:
            tp += 1
    for chrom in truth_dict:
        for pos in fp_dict[chrom]:
            fp += 1
    return tp, fp, fn



parser = argparse.ArgumentParser( description='Get metrics for xtea and tldr before and after runs')
parser.add_argument('--truth', required=True)
parser.add_argument('--tldr-before', required=True)
parser.add_argument('--xtea-before', required=True)
parser.add_argument('--somrit', required=True)
parser.add_argument('--window-size', type=int, default=500)
parser.add_argument('--fraction', type=float, required=True)
parser.add_argument('--rep', type=int, required=True)
parser.add_argument('--out-all', required=True)
parser.add_argument('--out-rt', required=True)
parser.add_argument('--chrom-list', default="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY")
args = parser.parse_args()

# Read in truth data for each repeat type
truth_all, truth_rt = get_truth_dict(args.truth, args.window_size, args.chrom_list.split(','))


tp_before, fp_before, fn_before = get_xtea_stats_all(args.xtea_before, truth_all)
if tp_before == 0 and fn_before == 0:
    xtea_before_pb_precision = 0
    xtea_before_pb_recall = 0
else:
    xtea_before_pb_precision = tp_before/(tp_before+len(fp_before))
    xtea_before_pb_recall = tp_before/(tp_before+fn_before)

out = str(args.fraction) + "\t" + str(args.rep) + "\t" + str(tp_before) + "\t" + str(len(fp_before)) + "\t" + str(fn_before) + "\t" + str(xtea_before_pb_precision) + "\t" + str(xtea_before_pb_recall) + "\t"

before_tp, before_fp, before_fn = get_tldr_stats_all(args.tldr_before, truth_all)
if before_tp == 0 and before_fn == 0:
    tldr_before_pb_precision = 0
    tldr_before_pb_recall = 0
else:
    tldr_before_pb_precision = before_tp/(before_tp+before_fp)
    tldr_before_pb_recall = before_tp/(before_tp+before_fn)

out += str(before_tp) + "\t" + str(before_fp) + "\t" + str(before_fn) + "\t" + str(tldr_before_pb_precision) + "\t" + str(tldr_before_pb_recall) + "\t"

tp, fp, fn = get_somrit_stats_all(args.somrit, args.window_size, truth_all)
precision = tp/(tp+fp)
recall = tp/(tp+fn)
out += str(tp) + "\t" + str(fp) + "\t" + str(fn) + "\t" + str(precision) + "\t" + str(recall)+"\n"
with open(args.out_all, 'w') as out_all:
    out_all.write(out)



tp_before, fp_before, fn_before = get_xtea_stats_rt(args.xtea_before, truth_rt)
if tp_before == 0 and fn_before == 0:
    xtea_before_pb_precision = 0
    xtea_before_pb_recall = 0
else:
    xtea_before_pb_precision = tp_before/(tp_before+len(fp_before))
    xtea_before_pb_recall = tp_before/(tp_before+fn_before)

out = str(args.fraction) + "\t" + str(args.rep) + "\t" + str(tp_before) + "\t" + str(len(fp_before)) + "\t" + str(fn_before) + "\t" + str(xtea_before_pb_precision) + "\t" + str(xtea_before_pb_recall) + "\t"

before_tp, before_fp, before_fn = get_tldr_stats_rt(args.tldr_before, truth_rt)
if before_tp == 0 and before_fn == 0:
    tldr_before_pb_precision = 0
    tldr_before_pb_recall = 0
else:
    tldr_before_pb_precision = before_tp/(before_tp+before_fp)
    tldr_before_pb_recall = before_tp/(before_tp+before_fn)

out += str(before_tp) + "\t" + str(before_fp) + "\t" + str(before_fn) + "\t" + str(tldr_before_pb_precision) + "\t" + str(tldr_before_pb_recall) + "\t"

tp, fp, fn = get_somrit_stats_rt(args.somrit, args.window_size, truth_rt)
if fp == 0 and tp == 0:
    precision = 0
else:
    precision = tp/(tp+fp)
if tp == 0 and fn == 0:
    recall = 0
else:
    recall = tp/(tp+fn)
out += str(tp) + "\t" + str(fp) + "\t" + str(fn) + "\t" + str(precision) + "\t" + str(recall)+"\n"
with open(args.out_rt, 'w') as out_rt:
    out_rt.write(out)

