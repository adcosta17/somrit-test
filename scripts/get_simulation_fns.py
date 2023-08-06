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

def get_locs(records_list):
    ret = ""
    for record in records_list:
        if ret != "":
            ret = ret + ","
        ret = ret + record.reference_name+":"+str(record.reference_start)+"-"+str(record.reference_end)
    return ret

def get_somrit_fns(file, window, truth_reads, bam_records):
    seen_reads = defaultdict(list)
    seen_fn = defaultdict(list)
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
            # Check if the position is within window size of a truth call
            match_name = match(row[3],truth_reads)
            if match_name is not None:
                if abs(truth_reads[match_name][2] - start) < window or abs(truth_reads[match_name][2] - end) < window:
                    seen_reads[match_name].append("\t".join(row)) 
                else:
                    seen_fn[match_name].append("\t".join(row))
    for read in truth_reads:
        if read in seen_reads:
            print("Found\t"+read+"\t"+"\t".join(list(map(str,truth_reads[read])))+"\t"+"\t".join(seen_reads[read]))
        elif read in seen_fn:
            alignment_locs = get_locs(bam_records[read])
            print("Missing_"+alignment_locs+"\t"+read+"\t"+"\t".join(list(map(str,truth_reads[read])))+"\t"+"\t".join(seen_fn[read]))
        else:
            if read in bam_records:
                # Have an alginment for the read that doesn't have an insert
                alignment_locs = get_locs(bam_records[read])
                print("Missing_"+alignment_locs+"\t"+read+"\t"+"\t".join(list(map(str,truth_reads[read]))))
            else:
                # Read not aligned at all
                print("Missing_No_Alignment\t"+read+"\t"+"\t".join(list(map(str,truth_reads[read]))))





parser = argparse.ArgumentParser( description='Get metrics for xtea and tldr before and after runs')
parser.add_argument('--inserts-tsv', required=True)
parser.add_argument('--spiked-reads', required=True)
parser.add_argument('--bam', required=True)
parser.add_argument('--somrit', required=True)
parser.add_argument('--window-size', type=int, default=1000)
parser.add_argument('--rep', type=int, required=True)
parser.add_argument('--chrom-list', default="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY")
args = parser.parse_args()

# Start by reading in the list of spiked in reads and the positions they come from
# Then read in the metadata on those positions that will tell us more about the positions
spike_in_reads = defaultdict(list)
spike_in_positions = defaultdict(IntervalTree)
insert_data = defaultdict(list)
bam_records = defaultdict(list)

with open(args.inserts_tsv, 'r') as in_tsv:
    for line in in_tsv:
        row = line.strip().split('\t')
        insert_data[row[0]] = row

samfile = pysam.AlignmentFile(args.bam)
for record in samfile.fetch():
    if record.query_name in insert_data:
        bam_records[record.query_name].append(record)

# Get the per read stats for somrit and tldr
# Positions for sniffles and xt
with open(args.spiked_reads, 'r') as in_reads:
    for line in in_reads:
        row = line.strip().split('\t')
        values = [row[0], insert_data[row[0]][5], int(insert_data[row[0]][6]), insert_data[row[0]][7], insert_data[row[0]][0], insert_data[row[0]][1], insert_data[row[0]][2], insert_data[row[0]][3], insert_data[row[0]][4]]
        spike_in_reads[row[1]] = values
        spike_in_positions[insert_data[row[0]][5]][int(insert_data[row[0]][6]):int(insert_data[row[0]][6])+1] = 1

# Get stats for each case
get_somrit_fns(args.somrit, args.window_size, spike_in_reads, bam_records)

