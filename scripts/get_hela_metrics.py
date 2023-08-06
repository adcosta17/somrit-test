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
import mappy as mp

def match(name, truth_reads):
    for read in truth_reads:
        if read in name:
            return read
    return None


def get_tldr_stats(file, line_start, line_end, mcherry_start, mcherry_end, mcherry_seq):
    seen = defaultdict(IntervalTree)
    count = 0
    a = mp.Aligner(mcherry_seq)  # load or build index
    if not a: raise Exception("ERROR: failed to load/build index")
    with open(file) as in_calls:
        for line in in_calls:
            if count == 0:
                count += 1
                continue
            row = line.strip().split('\t')
            if row[5] != "L1mCherry" or "NoTEAlignment" in row[24]:
                found = False
                for hit in a.map(row[21]): # traverse alignments
                    if hit.mapq >= 60:
                        if len(row[21]) >= 100 and (hit.r_en - hit.r_st) >= 100:
                            found = True
                    break
                if not found or int(row[15]) < 1 or int(row[16]) > 1:
                    continue
            else:
                if "LeftFlankSize" in row[24] or "RightFlankSize" in row[24] or "NoFamily" in row[24] or "NonRemappable" in row[24]:
                    continue
                if int(row[9]) < 100 or int(row[8])-int(row[7]) < 100 or int(row[15]) < 1 or int(row[16]) > 1:
                    continue
            # Check to see if the element overlaps the mcherry or maps outside the line
            if (int(row[7]) < mcherry_start and int(row[8]) > mcherry_start) or (int(row[7]) < mcherry_end and int(row[8]) > mcherry_end) or (int(row[7]) < line_end and int(row[8]) > line_end) or (int(row[7]) > line_end and int(row[8]) > mcherry_start):
                if int(row[2]) < int(row[3]):
                    seen[row[1]][int(row[2]):int(row[2])+1] = row
                else:
                    seen[row[1]][int(row[3]):int(row[3])+1] = row
    return seen


def get_somrit_stats(file, line_start, line_end, mcherry_start, mcherry_end, mcherry_seq):
    seen = defaultdict(IntervalTree)
    count = 0
    a = mp.Aligner(mcherry_seq)  # load or build index
    if not a: raise Exception("ERROR: failed to load/build index")
    with open(file) as in_calls:
        for line in in_calls:
            if count == 0:
                count += 1
                continue
            if "PASS" not in line:
                continue
            row = line.strip().split('\t')
            if "L1mCherry" not in row[7]:
                continue
            chrom = row[0]
            start = int(row[1])
            end = int(row[2])
            # Check how many samples we see
            items = row[5].split(',')
            samples = {}
            for item in items:
                sample = item.split(':')[0]
                samples[sample] = 1
            if len(samples) > 1:
                # See in more than one sample
                continue
            seen[chrom][start:end] = row
    final_counts = defaultdict(IntervalTree)
    count = 0
    with open(file) as in_calls:
        for line in in_calls:
            if count == 0:
                count += 1
                continue
            if "PASS" not in line or "Polymorphic" in line:
                continue
            row = line.strip().split('\t')
            if "L1mCherry" not in row[7]:
                continue
            if not ((int(row[14]) < mcherry_start and int(row[15]) > mcherry_start) or (int(row[14]) < mcherry_end and int(row[15]) > mcherry_end) or (int(row[14]) < line_end and int(row[15]) > line_end) or (int(row[14]) > line_end and int(row[15]) > mcherry_start)):
                found = False
                for hit in a.map(row[4]): # traverse alignments
                    if hit.mapq >= 60:
                        found = True
                    break
                if not found:
                    continue
            chrom = row[0]
            start = int(row[1])
            end = int(row[2])
            # Check how many samples we see
            items = row[5].split(',')
            samples = {}
            for item in items:
                sample = item.split(':')[0]
                samples[sample] = 1
            if len(samples) > 1:
                # See in more than one sample
                continue
            nearby = final_counts[chrom][start-2500:end+2500]
            if len(nearby) > 0:
                # Already added this position to the list
                continue
            nearby = seen[chrom][start-2500:end+2500]
            if len(nearby) == 1:
                final_counts[chrom][start:end] = [row]
            else:
                # have multiple insert supporting reads here, generate a list for the position
                ret = []
                samples = {}
                for item in nearby:
                    ret.append(item.data)
                    reads = item.data[5].split(',')
                    for r in reads:
                        sample = r.split(':')[0]
                        samples[sample] = 1
                if len(samples) > 1:
                    # Have calls that individually are unique to a sample but together are not
                    continue
                final_counts[chrom][start:end] = ret
    return final_counts



parser = argparse.ArgumentParser( description='Get metrics for xtea and tldr before and after runs')
parser.add_argument('--tldr', required=True)
parser.add_argument('--somrit', required=True)
parser.add_argument('--line-seq', required=True)
parser.add_argument('--l1mcherry-seq', required=True)
parser.add_argument('--mcherry-seq', required=True)
parser.add_argument('--insertions', required=True)
args = parser.parse_args()

count = 0
insert_positions = defaultdict(IntervalTree)
with open(args.insertions, 'r') as in_csv:
    for line in in_csv:
        if count == 0:
            count = 1
            continue
        row = line.strip().split(',')
        insert_positions[row[1]][int(row[2])-100:int(row[3])+100] = row[12].split('|')[0]

# Get the positions of the mcherry and LINE seqs in the L1mCherry seqeuence
a = mp.Aligner(args.l1mcherry_seq)  # load or build index
if not a: raise Exception("ERROR: failed to load/build index")
line_start = 0
line_end = 0
mcherry_start = 0
mcherry_end = 0
for name, seq, qual in mp.fastx_read(args.line_seq): # read a fasta/q sequence
    for hit in a.map(seq): # traverse alignments
        line_start = hit.r_st
        line_end = hit.r_en
for name, seq, qual in mp.fastx_read(args.mcherry_seq): # read a fasta/q sequence
    for hit in a.map(seq): # traverse alignments
        mcherry_start = hit.r_st
        mcherry_end = hit.r_en     

print(line_start)
print(line_end)
print(mcherry_start)
print(mcherry_end)

somrit_calls = get_somrit_stats(args.somrit, line_start, line_end, mcherry_start, mcherry_end, args.mcherry_seq)
tldr_calls = get_tldr_stats(args.tldr, line_start, line_end, mcherry_start, mcherry_end, args.mcherry_seq)


somrit_uniq = 0
tldr_uniq = 0
both = 0

in_set_tldr = 0
in_set_somrit = 0

somrit_single_read = 0

print("Chrom\tStart\tEnd\tSomirtSample\tSomritReadCount\tTldrSample\tTldrReadCount")
seen_tldr = defaultdict(IntervalTree)
for chrom in somrit_calls:
    for item in somrit_calls[chrom]:
        somrit_sample = item.data[0][5].split(':')[0]
        somrit_count = 0
        for row in item.data:
            reads = row[5].split(',')
            somrit_count += len(reads)
        if somrit_count == 1:
            somrit_single_read += 1
        nearby = tldr_calls[chrom][item.begin-2500:item.end+2500]
        tldr_sample = "NA"
        tldr_count = "NA"
        if len(nearby) > 0:
            # have a matching tldr call
            both += 1
            for r in nearby:
                seen_tldr[chrom][r.begin:r.end] = 1
                tldr_sample = r.data[17].split('|')[0]
                tldr_count = int(r.data[15])
                break
        else:
            somrit_uniq += 1
        nearby = insert_positions[chrom][item.begin:item.end]
        if len(nearby) > 0 :
            in_set_somrit += 1
        print(chrom+"\t"+str(item.begin)+"\t"+str(item.end)+"\t"+somrit_sample+"\t"+str(somrit_count)+"\t"+tldr_sample+"\t"+str(tldr_count))

# Now check tldr calls
for chrom in tldr_calls:
    for item in tldr_calls[chrom]:
        nearby = insert_positions[chrom][item.begin:item.end]
        if len(nearby) > 0 :
            in_set_tldr += 1
        nearby = seen_tldr[chrom][item.begin:item.end]
        if len(nearby) > 0:
            # Seen this with a somrit call
            continue
        # Call is unique to tldr
        tldr_sample = item.data[17].split('|')[0]
        tldr_count = int(item.data[15])
        print(chrom+"\t"+str(item.begin)+"\t"+str(item.end)+"\tNA\tNA\t"+tldr_sample+"\t"+str(tldr_count))
        tldr_uniq += 1

print("\n")
print(somrit_uniq)
print(tldr_uniq)
print(both)

print(in_set_somrit)
print(in_set_tldr)
print(somrit_single_read)
