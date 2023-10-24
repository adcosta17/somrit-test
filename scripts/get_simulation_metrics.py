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

def match(name, truth_reads):
    for read in truth_reads:
        if read in name:
            return read
    return None

def same_family(family, annotation):
    if "line" in family.lower() and "line" in annotation.lower():
        return True
    elif "alu" in family.lower() and "alu" in annotation.lower():
        return True
    elif "sva" in family.lower() and "sva" in annotation.lower():
        return True
    else:
        return False


def get_nearby(chrom, start, end, seen):
    if chrom in seen:
        count = 0
        i = start
        while i < end:
            if i in seen[chrom]:
                count += seen[chrom][i]
            i += 1
        if count >= 1:
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



def get_xtea_stats(file_prefix, window,truth_dict, min_size, max_size, family):
    seen = {}
    tp = 0
    fn = 0
    tp_distance = 0
    for chrom in truth_dict:
        seen[chrom] = defaultdict(int)
    with open(file_prefix+".merged_ALU.txt")as in_calls:
        for line in in_calls:
            row = line.strip().split('\t')
            chrom = row[0]
            start = int(row[1])
            seq_len = len(row[9])
            if family != "all":
                if not same_family(family, row[2]):
                    continue
            if seq_len < min_size or seq_len > max_size:
                continue
            # Check if the position is within window size of a truth call
            nearby = truth_dict[chrom][start-window:start+window]
            if len(nearby) > 0:
                # Have a tp
                for item in nearby:
                    #print(item.data)
                    if same_family(row[2], item.data[0]):
                        if int(item.begin) in seen[chrom]:
                            seen[chrom][int(item.begin)] = min(abs(int(item.begin)-start),seen[chrom][int(item.begin)])
                        else:
                            seen[chrom][int(item.begin)] = abs(int(item.begin)-start)
    with open(file_prefix+".merged_LINE1.txt")as in_calls:
        for line in in_calls:
            row = line.strip().split('\t')
            chrom = row[0]
            start = int(row[1])
            seq_len = len(row[9])
            if family != "all":
                if not same_family(family, row[2]):
                    continue
            if seq_len < min_size or seq_len > max_size:
                continue
            # Check if the position is within window size of a truth call
            nearby = truth_dict[chrom][start-window:start+window]
            if len(nearby) > 0:
                # Have a tp
                for item in nearby:
                    #print(item.data)
                    if same_family(row[2], item.data[0]):
                        if int(item.begin) in seen[chrom]:
                            seen[chrom][int(item.begin)] = min(abs(int(item.begin)-start),seen[chrom][int(item.begin)])
                        else:
                            seen[chrom][int(item.begin)] = abs(int(item.begin)-start)
    with open(file_prefix+".merged_HERV.txt")as in_calls:
        for line in in_calls:
            row = line.strip().split('\t')
            chrom = row[0]
            start = int(row[1])
            seq_len = len(row[9])
            if family != "all":
                if not same_family(family, row[2]):
                    continue
            if seq_len < min_size or seq_len > max_size:
                continue
            # Check if the position is within window size of a truth call
            nearby = truth_dict[chrom][start-window:start+window]
            if len(nearby) > 0:
                # Have a tp
                for item in nearby:
                    #print(item.data)
                    if same_family(row[2], item.data[0]):
                        if int(item.begin) in seen[chrom]:
                            seen[chrom][int(item.begin)] = min(abs(int(item.begin)-start),seen[chrom][int(item.begin)])
                        else:
                            seen[chrom][int(item.begin)] = abs(int(item.begin)-start)
    with open(file_prefix+".merged_SVA.txt")as in_calls:
        for line in in_calls:
            row = line.strip().split('\t')
            chrom = row[0]
            start = int(row[1])
            seq_len = len(row[9])
            if family != "all":
                if not same_family(family, row[2]):
                    continue
            if seq_len < min_size or seq_len > max_size:
                continue
            # Check if the position is within window size of a truth call
            nearby = truth_dict[chrom][start-window:start+window]
            if len(nearby) > 0:
                # Have a tp
                for item in nearby:
                    #print(item.data)
                    if same_family(row[2], item.data[0]):
                        if int(item.begin) in seen[chrom]:
                            seen[chrom][int(item.begin)] = min(abs(int(item.begin)-start),seen[chrom][int(item.begin)])
                        else:
                            seen[chrom][int(item.begin)] = abs(int(item.begin)-start)
    for chrom in truth_dict:
        for item in truth_dict[chrom]:
            if int(item.begin) not in seen[chrom]:
                fn += 1
    for chrom in seen:
        for pos in seen[chrom]:
            tp += 1
            tp_distance += seen[chrom][pos]
    return tp, fn


def get_tldr_stats(file, window, truth_reads, truth_dict, min_size, max_size, family):
    seen = {}
    seen_reads = {}
    tp = 0
    fn = 0
    tp_reads = 0
    fn_reads = 0
    fp = 0
    for chrom in truth_dict:
        seen[chrom] = defaultdict(int)
    count = 0
    with open(file) as in_calls:
        for line in in_calls:
            if count == 0:
                count += 1
                continue
            row = line.strip().split('\t')
            if "PASS" not in line:
                continue
            chrom = row[1]
            start = int(row[2])
            end = int(row[3])
            seq_len = int(row[9])
            if family != "all":
                if not same_family(family, row[5]):
                    continue
            if seq_len < min_size or seq_len > max_size:
                continue
            # Check if the position is within window size of a truth call
            nearby = truth_dict[chrom][start-window:end+window]
            if len(nearby) > 0:
                # Have a tp
                found = False
                for item in nearby:
                    #print(item.data)
                    if same_family(item.data[0], row[5]):
                        seen[chrom][int(item.begin)] = 1
                        if row[0] in truth_reads:
                            if abs(truth_reads[row[0]][2] - start) < window or abs(truth_reads[row[0]][2] - end) < window:
                                seen_reads[row[0]] = 1
                                found = True
                if not found:
                    fp += 1
            else:
                fp += 1
    for chrom in truth_dict:
        for item in truth_dict[chrom]:
            if int(item.begin) not in seen[chrom]:
                fn += 1
    for read in truth_reads:
        if read not in seen_reads:
            fn_reads += 1
    for chrom in seen:
        for pos in seen[chrom]:
            tp += 1
    tp_reads = len(seen_reads)
    return tp, fn, tp_reads, fn_reads, fp


def get_somrit_stats(file, window, truth_reads, truth_dict, min_size, max_size, family):
    base_tree = defaultdict(IntervalTree)
    fp_dict = defaultdict(IntervalTree)
    seen = {}
    all_inserts = {}
    tp = 0
    fn = 0
    fp = 0
    fn_reads = 0
    found_count = defaultdict(int)
    seen_reads = {}
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
            seq_len = len(row[4])
            if family != "all":
                if not same_family(family, row[7]):
                    continue
            if seq_len < min_size or seq_len > max_size:
                continue
            # Check if the position is within window size of a truth call
            nearby = truth_dict[chrom][start-window:end+window]
            if len(nearby) > 0:
                # Have a tp
                for item in nearby:
                    #print(item.data)
                    if same_family(item.data[0], row[7]):
                        seen[chrom][int(item.begin)] = 1
                        match_name = match(row[3],truth_reads)
                        if match_name is not None:
                            if abs(truth_reads[match_name][2] - start) < window or abs(truth_reads[match_name][2] - end) < window:
                                seen_reads[match_name] = 1 
            if chrom not in all_inserts:
                all_inserts[chrom] = defaultdict(int)
            all_inserts[chrom][start] = 1
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
            if family != "all":
                if not same_family(family, row[7]):
                    continue
            if seq_len < min_size or seq_len > max_size:
                continue
            # Check if the position is within window size of a truth call
            nearby = truth_dict[chrom][start-window:end+window]
            found = False
            if len(nearby) > 0:
                # Have a tp
                for item in nearby:
                    #print(item.data)
                    if same_family(item.data[0], row[7]):
                        found = True
            if found:
                continue
            nearby = get_nearby(chrom, start-window, end+window, all_inserts)
            nearby_fp = fp_dict[chrom][start:end]
            if nearby and len(nearby_fp) == 0:
                fp_dict[chrom][(start-window):(end+window)] = start
    for chrom in truth_dict:
        for item in truth_dict[chrom]:
            if int(item.begin) not in seen[chrom]:
                fn += 1
    for read in truth_reads:
        if read not in seen_reads:
            fn_reads += 1
    for chrom in seen:
        for pos in seen[chrom]:
            tp += 1
    tp_reads = len(seen_reads)
    for chrom in all_inserts:
        for pos in fp_dict[chrom]:
            fp += 1
    return tp, fn, tp_reads, fn_reads, fp


def get_somrit_fns(file, window, truth_reads, min_size, max_size, family):
    seen_reads = {}
    seen_fn = {}
    tp_reads = 0
    fn_reads = 0
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
            seq_len = len(row[4])
            if family != "all":
                if not same_family(family, row[7]):
                    continue
            if seq_len < min_size or seq_len > max_size:
                continue
            # Check if the position is within window size of a truth call
            match_name = match(row[3],truth_reads)
            if match_name is not None:
                if abs(truth_reads[match_name][2] - start) < window or abs(truth_reads[match_name][2] - end) < window:
                    seen_reads[match_name] = row 
                else:
                    seen_fn[match_name] = row
    for read in truth_reads:
        if read in seen_reads:
            print("Found\t"+read+"\t"+"\t".join(list(map(str,truth_reads[read])))+"\t"+"\t".join(seen_reads[read]))
        elif read in seen_fn:
            print("Missing\t"+read+"\t"+"\t".join(list(map(str,truth_reads[read])))+"\t"+"\t".join(seen_fn[read]))
        else:
            print("Missing\t"+read+"\t"+"\t".join(list(map(str,truth_reads[read]))))

def get_sniffles_stats(file, window, truth_dict, min_size, max_size, family):
    seen = {}
    tp = 0
    fn = 0
    for chrom in truth_dict:
        seen[chrom] = defaultdict(int)
    bcf_in = pysam.VariantFile(file)
    for rec in bcf_in.fetch():
        if rec.info["SVTYPE"] == "INS":
            chrom = rec.contig
            start = rec.pos
            seq_len = len(rec.alts[0])
            if family != "all":
                if not same_family(family, "NA"):
                    continue
            if seq_len < min_size or seq_len > max_size:
                continue
            nearby = truth_dict[chrom][start-window:start+window]
            if len(nearby) > 0:
                # Have a tp
                for item in nearby:
                    #print(item.data)
                    if int(item.begin) in seen[chrom]:
                        seen[chrom][int(item.begin)] = min(abs(int(item.begin)-start),seen[chrom][int(item.begin)])
                    else:
                        seen[chrom][int(item.begin)] = abs(int(item.begin)-start)
    for chrom in truth_dict:
        for item in truth_dict[chrom]:
            if int(item.begin) not in seen[chrom]:
                fn += 1
    for chrom in seen:
        for pos in seen[chrom]:
            tp += 1
    return tp, fn






parser = argparse.ArgumentParser( description='Get metrics for xtea and tldr before and after runs')
parser.add_argument('--inserts-tsv', required=True)
parser.add_argument('--spiked-reads', required=True)
parser.add_argument('--tldr', required=True)
parser.add_argument('--tldr-realign', required=True)
parser.add_argument('--xtea', required=True)
parser.add_argument('--somrit', required=True)
parser.add_argument('--sniffles', required=True)
parser.add_argument('--min-size', default=0, type=int)
parser.add_argument('--max-size', default=10000, type=int)
parser.add_argument('--family', default="all")
parser.add_argument('--window-size', type=int, default=500)
parser.add_argument('--rep', type=int, required=True)
parser.add_argument('--chrom-list', default="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY")
args = parser.parse_args()

# Start by reading in the list of spiked in reads and the positions they come from
# Then read in the metadata on those positions that will tell us more about the positions
spike_in_reads = defaultdict(list)
spike_in_positions = defaultdict(IntervalTree)
insert_data = defaultdict(list)

with open(args.inserts_tsv, 'r') as in_tsv:
    for line in in_tsv:
        row = line.strip().split('\t')
        insert_data[row[0]] = row

# Get the per read stats for somrit and tldr
# Positions for sniffles and xt
with open(args.spiked_reads, 'r') as in_reads:
    for line in in_reads:
        row = line.strip().split('\t')
        values = [row[0], insert_data[row[0]][5], int(insert_data[row[0]][6]), insert_data[row[0]][7], insert_data[row[0]][0], insert_data[row[0]][1], insert_data[row[0]][2], insert_data[row[0]][3], insert_data[row[0]][4]]
        if len(insert_data[row[0]][8]) < args.min_size or len(insert_data[row[0]][8]) > args.max_size:
            continue
        if args.family != "all":
            if not same_family(args.family, insert_data[row[0]][7]):
                continue
        spike_in_reads[row[1]] = values
        spike_in_positions[insert_data[row[0]][5]][int(insert_data[row[0]][6]):int(insert_data[row[0]][6])+1] = [insert_data[row[0]][7], len(insert_data[row[0]][8])]

# Get stats for each case
tp_somrit, fn_somrit, tp_reads_somrit, fn_reads_somrit, fp_somrit = get_somrit_stats(args.somrit,args.window_size,spike_in_reads,spike_in_positions, args.min_size, args.max_size, args.family)
tp_tldr, fn_tldr, tp_reads_tldr, fn_reads_tldr, fp_tldr = get_tldr_stats(args.tldr,args.window_size,spike_in_reads,spike_in_positions, args.min_size, args.max_size, args.family)
tp_tldr_realign, fn_tldr_realign, tp_reads_tldr_realign, fn_reads_tldr_realign, fp_tldr_realign = get_tldr_stats(args.tldr_realign,args.window_size,spike_in_reads,spike_in_positions, args.min_size, args.max_size, args.family)
tp_xtea, fn_xtea = get_xtea_stats(args.xtea,args.window_size,spike_in_positions, args.min_size, args.max_size, args.family)
tp_sniffles, fn_sniffiles = get_sniffles_stats(args.sniffles,args.window_size,spike_in_positions, args.min_size, args.max_size, args.family)
if args.rep > 0:
    to_print = [str(args.rep), str(tp_somrit), str(fn_somrit), str(float(tp_somrit)/(tp_somrit+fn_somrit)), str(tp_reads_somrit), str(fn_reads_somrit), str(float(tp_reads_somrit)/(tp_reads_somrit+fn_reads_somrit)), str(tp_tldr), str(fn_tldr), str(float(tp_tldr)/(tp_tldr+fn_tldr)), str(tp_reads_tldr), str(fn_reads_tldr), str(float(tp_reads_tldr)/(tp_reads_tldr+fn_reads_tldr)), str(tp_xtea), str(fn_xtea), str(float(tp_xtea)/(tp_xtea+fn_xtea)), str(tp_sniffles), str(fn_sniffiles), str(float(tp_sniffles)/(tp_sniffles+fn_sniffiles)),str(tp_tldr_realign), str(fn_tldr_realign), str(float(tp_tldr_realign)/(tp_tldr_realign+fn_tldr_realign)), str(fp_somrit), str(float(tp_somrit/(tp_somrit+fp_somrit))), str(fp_tldr), str(float(tp_tldr/(tp_tldr+fp_tldr))), str(fp_tldr_realign), str(float(tp_tldr_realign/(tp_tldr_realign+fp_tldr_realign)))]
else:
    to_print = [str(args.rep), str(fp_somrit), str(fp_tldr)]
print("\t".join(to_print))
#get_somrit_fns(args.somrit,args.window_size,spike_in_reads,args.min_size, args.max_size, args.family)

