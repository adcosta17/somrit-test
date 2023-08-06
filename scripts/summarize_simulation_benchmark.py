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



parser = argparse.ArgumentParser( description='Get metrics for xtea and tldr before and after runs')
parser.add_argument('--rep', type=int, required=True)
parser.add_argument('--sample', type=int, required=True)
parser.add_argument('--chroms', default="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY")
args = parser.parse_args()


sniffles_before_mem = 0
sniffles_before_time = 0
with open("benchmark/sniffles/HG"+str(args.sample)+".spike_in."+str(args.rep)+".txt", 'r') as in_file:
    count = 0
    for line in in_file:
        if count == 0:
            count = 1
            continue
        row = line.split("\t")
        sniffles_before_time = row[0]
        sniffles_before_mem = row[2]


xtea_before_mem = 0
xtea_before_time = 0
with open("benchmark/xtea/HG"+str(args.sample)+".spike_in."+str(args.rep)+".txt", 'r') as in_file:
    count = 0
    for line in in_file:
        if count == 0:
            count = 1
            continue
        row = line.split("\t")
        xtea_before_time = row[0]
        xtea_before_mem = row[2]


tldr_before_mem = 0
tldr_before_time = 0
with open("benchmark/tldr/HG"+str(args.sample)+".spike_in."+str(args.rep)+".txt", 'r') as in_file:
    count = 0
    for line in in_file:
        if count == 0:
            count = 1
            continue
        row = line.split("\t")
        tldr_before_time = row[0]
        tldr_before_mem = row[2]


somrit_extract_mem = 0
somrit_extract_time = 0
with open("benchmark/extract/HG"+str(args.sample)+".spike_in."+str(args.rep)+".txt", 'r') as in_file:
    count = 0
    for line in in_file:
        if count == 0:
            count = 1
            continue
        row = line.split("\t")
        somrit_extract_time = row[0]
        somrit_extract_mem = row[2]

somrit_realign_mem = 0
somrit_realign_time = 0
somrit_realign_mem_max = 0
somrit_realign_time_max = 0
for chrom in args.chroms.split(','):
    with open("benchmark/realign/HG"+str(args.sample)+"."+chrom+".spike_in."+str(args.rep)+".txt", 'r') as in_file:
        count = 0
        for line in in_file:
            if count == 0:
                count = 1
                continue
            row = line.split("\t")
            somrit_realign_time += float(row[0])
            somrit_realign_mem += float(row[2])
            if float(row[0]) > somrit_realign_time_max:
                somrit_realign_time_max = float(row[0])
            if float(row[2]) > somrit_realign_mem_max:
                somrit_realign_mem_max = float(row[2])

somrit_realign_time_avg = somrit_realign_time/len(args.chroms.split(','))
somrit_realign_mem_avg = somrit_realign_mem/len(args.chroms.split(','))


somrit_merge_mem = 0
somrit_merge_time = 0
with open("benchmark/merge/HG"+str(args.sample)+".spike_in."+str(args.rep)+".txt", 'r') as in_file:
    count = 0
    for line in in_file:
        if count == 0:
            count = 1
            continue
        row = line.split("\t")
        somrit_merge_time = row[0]
        somrit_merge_mem = row[2]


somrit_classified_mem = 0
somrit_classified_time = 0
with open("benchmark/classify/HG"+str(args.sample)+".spike_in."+str(args.rep)+".txt", 'r') as in_file:
    count = 0
    for line in in_file:
        if count == 0:
            count = 1
            continue
        row = line.split("\t")
        somrit_classified_time = row[0]
        somrit_classified_mem = row[2]


somrit_filtered_mem = 0
somrit_filtered_time = 0
with open("benchmark/filter/HG"+str(args.sample)+".spike_in."+str(args.rep)+".txt", 'r') as in_file:
    count = 0
    for line in in_file:
        if count == 0:
            count = 1
            continue
        row = line.split("\t")
        somrit_filtered_time = row[0]
        somrit_filtered_mem = row[2]

# Print out values
print(str(args.sample)+"\t0\t"+str(args.rep)+"\t0\t0\t"+str(sniffles_before_mem)+"\t"+str(sniffles_before_time)+"\t"+str(xtea_before_mem)+"\t"+str(xtea_before_time)+"\t"+str(tldr_before_mem)+"\t"+str(tldr_before_time)+"\t"+str(somrit_extract_mem)+"\t"+str(somrit_extract_time)+"\t"+str(somrit_realign_mem)+"\t"+str(somrit_realign_time)+"\t"+str(somrit_realign_mem_max)+"\t"+str(somrit_realign_time_max)+"\t"+str(somrit_realign_mem_avg)+"\t"+str(somrit_realign_time_avg)+"\t"+str(somrit_merge_mem)+"\t"+str(somrit_merge_time)+"\t"+str(somrit_classified_mem)+"\t"+str(somrit_classified_time)+"\t"+str(somrit_filtered_mem)+"\t"+str(somrit_filtered_time))
