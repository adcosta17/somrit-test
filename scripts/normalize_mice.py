import pysam
import argparse
import sys
import csv
import os
import math
from os import listdir
from os.path import isfile, join
from intervaltree import Interval, IntervalTree

from collections import defaultdict

def get_nearby(chrom, start, end, seen):
    i = start
    count = 0
    nearby = {}
    if chrom not in seen:
        return nearby
    while i < end:
        if i in seen[chrom]:
            nearby[i] = seen[chrom][i]
            count += 1
        i += 1
    return nearby

def get_somrit_test(file, test_sample, parental_calls):
    seen = {}
    insert_seen = {}
    count = 0
    with open(file) as in_calls:
        for line in in_calls:
            if count == 0:
                count += 1
                continue
            count += 1
            if count % 100000 == 0:
                print(count,file=sys.stderr)
            row = line.strip().split('\t')
            if "PASS" not in line or (int(row[13])-int(row[12]))/int(len(row[4])) < 0.5:
                continue
            chrom = row[0]
            start = int(row[1])
            end = int(row[2])
            if chrom not in seen:
                seen[chrom] = defaultdict(list)
            i = start
            while i < end:
                seen[chrom][i].append(row)
                i += 1
    final_counts = defaultdict(IntervalTree)
    count = 0
    with open(file) as in_calls:
        for line in in_calls:
            if count == 0:
                count += 1
                continue
            count += 1
            if count % 100000 == 0:
                print(count,file=sys.stderr)
            row = line.strip().split('\t')
            if "PASS" not in line or "Polymorphic" in line or "No_TSD_Found" in line or "Has_Poly_A_Tail" not in line or (int(row[13])-int(row[12]))/int(len(row[4])) < 0.5:
                continue
            chrom = row[0]
            start = int(row[1])
            end = int(row[2])
            final_nearby = final_counts[chrom][start-1000:end+1000]
            parental_nearby = parental_calls[chrom][start:end]
            if len(parental_nearby) > 0 or len(final_nearby) > 0:
                continue
            # Check how many samples we see
            if test_sample not in row[5]:
                continue
            nearby = get_nearby(chrom, start-1000, end+1000, seen)
            if len(nearby) == 1:
                samples = {}
                reads = row[5].split(',')
                for r in reads:
                    sample = r.split(':')[0]
                    samples[sample] = 1
                if len(samples) > 1 :
                    # Have calls that individually are unique to a sample but together are not
                    continue
                final_counts[chrom][start:end] = [row]
                #print(samples)
                #print(control_samples)
                #print([row])
            else:
                # have multiple insert supporting reads here, generate a list for the position
                ret = []
                samples = {}
                for i in nearby:
                    for item in nearby[i]:
                        if "PASS" in item[6]:
                            ret.append(item)
                        reads = item[5].split(',')
                        for r in reads:
                            sample = r.split(':')[0]
                            samples[sample] = 1
                if len(samples) > 1:
                    # Have calls that individually are unique to a sample but together are not
                    continue
                #print(samples)
                if len(ret) > 0:
                    final_counts[chrom][start:end] = ret
                #print(ret)
    return final_counts


parser = argparse.ArgumentParser( description='Take in a summary counts file for a sample and normalize it using effective bases')
parser.add_argument('--sample', required=True)
parser.add_argument('--bam', required=True)
parser.add_argument('--tsv', required=True)
parser.add_argument('--truth-csv', required=True)
parser.add_argument('--min-mapq', default=20, type=int)
parser.add_argument('--flank', default=1000, type=int)
parser.add_argument('--parental', required=True)
args = parser.parse_args()

count = 0
parental_calls = defaultdict(IntervalTree)
with open(args.parental, 'r') as in_parental:
    for line in in_parental:
        if count == 0:
            count = 1
            continue
        row = line.strip().split('\t')
        if len(row) > 4:
            parental_calls[row[2]][int(row[4])-1000:int(row[3])+1000] = 1

reads_to_use = {}
effective_bases_line = 0
effective_bases_alu = 0
effective_bases_sva = 0
effective_bases_erv = 0
effective_bases_all = 0
# Open bam and get a list of reads that map with at least mapq >= minimum_mapping_quality
count = 0
header = pysam.AlignmentFile(args.bam).header
for sq in header['SQ']:
    #print(sq['SN'])
    sam_reader = pysam.AlignmentFile(args.bam)
    tmp_sam_reader = sam_reader.fetch(contig=sq['SN'])
    for record in tmp_sam_reader:
        count += 1
        if count % 100000 == 0:
            print(count,file=sys.stderr)
        if record.mapping_quality >= args.min_mapq:
            if record.query_name not in reads_to_use and not record.is_secondary and not record.is_supplementary:
                reads_to_use[record.query_name] = record.query_alignment_length

line_sizes = []
alu_sizes = []
erv_sizes = []
sva_sizes = []
all_sizes = []
# Compute averaage insert sizes
with open(args.tsv) as csvfile:
    count = 0
    for row in csvfile:
        row_args = row.strip().split("\t")
        if count == 0:
            count = 1
            continue
        if "PASS" in row_args[8] and "ambiguous" not in row:
            count += 1
            all_sizes.append(len(row_args[7]))
            if "LINE" in row:
                line_sizes.append(len(row_args[7]))
            if "SINE" in row:
                alu_sizes.append(len(row_args[7]))
            if "SVA" in row:
                sva_sizes.append(len(row_args[7]))
            if "ERV" in row:
                erv_sizes.append(len(row_args[7]))

#print("Count: "+str(count))
avg_line_size = 0
if len(line_sizes) > 0:
    avg_line_size = sum(line_sizes)/len(line_sizes)
avg_alu_size = 0
if len(alu_sizes) > 0:
    avg_alu_size = sum(alu_sizes)/len(alu_sizes)
avg_erv_size = 0
if len(erv_sizes) > 0:
    avg_erv_size = sum(erv_sizes)/len(erv_sizes)
avg_sva_size = 0
if len(sva_sizes) > 0:
    avg_sva_size = sum(sva_sizes)/len(sva_sizes)
avg_all_size = 0
if len(all_sizes) > 0:
    avg_all_size = sum(all_sizes)/len(all_sizes)

#print(avg_line_size)
#print(avg_alu_size)
#print(avg_sva_size)
#print(avg_erv_size)

for read in reads_to_use:
    e_line = reads_to_use[read] - (2*args.flank + avg_line_size)
    if e_line > 0:
        effective_bases_line += e_line
    e_alu = reads_to_use[read] - (2*args.flank + avg_alu_size)
    if e_alu > 0:
        effective_bases_alu += e_alu
    e_sva = reads_to_use[read] - (2*args.flank + avg_sva_size)
    if e_sva > 0:
        effective_bases_sva += e_sva
    e_erv = reads_to_use[read] - (2*args.flank + avg_erv_size)
    if e_erv > 0:
        effective_bases_erv += e_erv
    e_all = reads_to_use[read] - (2*args.flank + avg_all_size)
    if e_all > 0:
        effective_bases_all += e_all

effective_bases_line = effective_bases_line/1000000000
effective_bases_alu = effective_bases_alu/1000000000
effective_bases_sva = effective_bases_sva/1000000000
effective_bases_erv = effective_bases_erv/1000000000
effective_bases_all = effective_bases_all/1000000000

#print(effective_bases_line)
#print(effective_bases_alu)
#print(effective_bases_sva)
#print(effective_bases_erv)


passing_line = 0
passing_alu = 0
passing_erv = 0
passing_sva = 0
total = 0
somrit_calls = get_somrit_test(args.tsv, args.sample, parental_calls)
print(somrit_calls, file=sys.stderr)
for chrom in somrit_calls:
    for item in somrit_calls[chrom]:
        line_count = 0
        alu_count = 0
        sva_count = 0
        erv_count = 0
        total += 1
        for row in item.data:
            if "l1" in row[7].lower():
                line_count += 1
            if "sine" in row[7].lower():
                alu_count += 1
            if "sva" in row[7].lower():
                sva_count += 1
            if "erv" in row[7].lower():
                erv_count += 1
        if line_count > 0 and alu_count == 0 and sva_count == 0 and erv_count == 0:
            passing_line += 1
        elif line_count == 0 and alu_count > 0 and sva_count == 0 and erv_count == 0:
            passing_alu += 1
        elif line_count == 0 and alu_count == 0 and sva_count > 0 and erv_count == 0:
            passing_sva += 1
        elif line_count == 0 and alu_count == 0 and sva_count == 0 and erv_count > 0:
            passing_erv += 1



print("Sample\tLINE\tLINE_Effective\tALU\tALU_Effective\tSVA\tSVA_Effective\tERV_Effective\tERV\tERV_Novel")
out = args.sample+"\t"
out += str(total)+"\t"+str(total/effective_bases_all)+"\t"
out += str(passing_line)+"\t"+str(passing_line/effective_bases_line)+"\t"
out += str(passing_alu)+"\t"+str(passing_alu/effective_bases_alu)+"\t"
out += str(passing_sva)+"\t"+str(passing_sva/effective_bases_sva)+"\t"
out += str(passing_erv)+"\t"+str(passing_erv/effective_bases_erv)
print(out)

count = 0
sample_count = 0
with open(args.truth_csv , 'r') as in_csv:
    for line in in_csv:
        if count == 0:
            count = 1
            continue
        row = line.strip().split(',')
        if len(row) < 16:
            continue
        if args.sample in row[15]:
            sample_count += 1

print(args.sample+"\t"+str(sample_count)+"\t"+str(sample_count/effective_bases_all))

