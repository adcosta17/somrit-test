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

def get_read_distance(pos, record):
    ref_count = record.reference_start
    read_count = 0
    in_window = False
    min_distance = 500
    for cg in re.findall('[0-9]*[A-Z=]', record.cigarstring):
        if cg.endswith('M'):
            ref_count += int(cg[:cg.find("M")])
            read_count += int(cg[:cg.find("M")])
        elif cg.endswith('X'):
            ref_count += int(cg[:cg.find("X")])
            read_count += int(cg[:cg.find("X")])
        elif cg.endswith('='):
            ref_count += int(cg[:cg.find("=")])
            read_count += int(cg[:cg.find("=")])
        elif cg.endswith('I'):
            read_count += int(cg[:cg.find("I")])
            if int(cg[:cg.find("I")]) > 50 and in_window:
                min_distance = min(min_distance, abs(ref_count-pos))
        elif cg.endswith('S'):
            read_count += int(cg[:cg.find("S")])
        elif cg.endswith('H'):
            read_count += int(cg[:cg.find("H")])
        elif cg.endswith('D'):
            ref_count += int(cg[:cg.find("D")])
        if abs(ref_count - pos) <= 500:
            in_window = True
        else:
            in_window = False
    return min_distance

def get_supporting_insertion_positions(chrom, start, bam_in, read_list):
    distances = {}
    for record in bam_in.fetch(contig=chrom, start=start-500, end=start+500):
        if read_list is not None:
            if record.query_name not in read_list: 
                continue
        distances[record.query_name] = str(get_read_distance(start, record))
    return distances

def get_vcf_distance(file, truth_dict, min_size, bam):
    seen = {}
    for chrom in truth_dict:
        seen[chrom] = {}
    bcf_in = pysam.VariantFile(file)
    bam_in = pysam.AlignmentFile(bam)
    for rec in bcf_in.fetch():
        if rec.info["SVTYPE"] == "INS":
            if len(rec.alts[0]) < min_size:
                continue
            chrom = rec.contig
            if chrom != "chr2":
                continue
            start = rec.pos
            nearby = truth_dict[chrom][start:start+1]
            closest = 999999999
            pos = -1
            if len(nearby) > 0:
                # If we actually have a nearby truth position, get the list of reads that support this SV and the nearest insertion position or softclip to expected position and the distance to it
                read_list = list(rec.info["RNAMES"])
                # Get the closest truth pos for this insert, make sure it only counts towards one
                for item in nearby:
                    if abs(int(item.data)-start) < closest:
                        pos = int(item.data)
                if pos > -1:
                    #print(item.data)
                    if pos in seen[chrom]:
                        if abs(pos-start) < seen[chrom][pos][0]:
                            read_distances = get_supporting_insertion_positions(chrom, start, bam_in, read_list)
                            seen[chrom][pos] = [abs(pos-start), read_distances, len(read_list)]
                    else:
                        read_distances = get_supporting_insertion_positions(chrom, start, bam_in, read_list)
                        seen[chrom][pos] = [abs(pos-start), read_distances, len(read_list)]
    return seen

def get_distance_diff(seen_before, seen_after):
    ret_before = []
    ret_after = []
    for chrom in seen_before:
        for pos in seen_before[chrom]:
            if pos in seen_after[chrom]:
                # Check each read
                for read in seen_before[chrom][pos][1]:
                    if read in seen_after[chrom][pos][1]:
                        ret_before.append(seen_before[chrom][pos][1][read])
                        ret_after.append(seen_after[chrom][pos][1][read])
    return ",".join(ret_before), ",".join(ret_after)
                
def get_support_diff(seen_before, seen_after):
    ret_before = []
    ret_after = []
    for chrom in seen_before:
        for pos in seen_before[chrom]:
            if pos in seen_after[chrom]:
                ret_before.append(str(seen_before[chrom][pos][2]))
                ret_after.append(str(seen_after[chrom][pos][2]))
    return ",".join(ret_before), ",".join(ret_after)


parser = argparse.ArgumentParser( description='Get metrics for Sniffles and CuteSV before and after runs')
parser.add_argument('--truth-mat', required=True)
parser.add_argument('--truth-pat', required=True)
parser.add_argument('--sniffles-before', required=True)
parser.add_argument('--sniffles-after', required=True)
parser.add_argument('--cutesv-before', required=True)
parser.add_argument('--cutesv-after', required=True)
parser.add_argument('--base-bam', required=True)
parser.add_argument('--realigned-bam', required=True)
parser.add_argument('--window-size', type=int, default=1000)
parser.add_argument('--fraction', type=float, required=True)
parser.add_argument('--rep', type=int, required=True)
parser.add_argument('--chrom-list', default="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY")
parser.add_argument('--insert-stats', required=True)
parser.add_argument('--rt-insert-stats', required=True)
args = parser.parse_args()

# Read in truth data for each repeat type
truth_large = get_truth_dict(args.truth_mat, args.truth_pat, args.window_size, args.chrom_list.split(','), 500)
truth_base = get_truth_dict(args.truth_mat, args.truth_pat, args.window_size, args.chrom_list.split(','), 50)

with open(args.insert_stats, 'w') as out_is, open(args.rt_insert_stats, 'w') as out_rt:
    seen_before = get_vcf_distance(args.sniffles_before, truth_base, 50, args.base_bam)
    seen_after = get_vcf_distance(args.sniffles_after, truth_base, 50, args.realigned_bam)
    diff_before_sniffles_50, diff_after_sniffles_50 = get_distance_diff(seen_before, seen_after)
    support_before_sniffles_50, support_after_sniffles_50 = get_support_diff(seen_before, seen_after)
    out_is.write(str(args.fraction)+"\t"+str(args.rep)+"\t"+support_before_sniffles_50+"\t"+support_after_sniffles_50+"\t"+diff_before_sniffles_50+"\t"+diff_after_sniffles_50+"\t")

    seen_before = get_vcf_distance(args.sniffles_before, truth_large, 500, args.base_bam)
    seen_after = get_vcf_distance(args.sniffles_after, truth_large, 500, args.realigned_bam)
    diff_before_sniffles_500, diff_after_sniffles_500 = get_distance_diff(seen_before, seen_after)
    support_before_sniffles_500, support_after_sniffles_500 = get_support_diff(seen_before, seen_after)
    out_rt.write(str(args.fraction)+"\t"+str(args.rep)+"\t"+support_before_sniffles_500+"\t"+support_after_sniffles_500+"\t"+diff_before_sniffles_500+"\t"+diff_after_sniffles_500+"\t")

    seen_before = get_vcf_distance(args.cutesv_before, truth_base, 50, args.base_bam)
    seen_after = get_vcf_distance(args.cutesv_after, truth_base, 50, args.realigned_bam)
    diff_before_cutesv_50, diff_after_cutesv_50 = get_distance_diff(seen_before, seen_after)
    support_before_cutesv_50, support_after_cutesv_50 = get_support_diff(seen_before, seen_after)
    out_is.write(support_before_cutesv_50+"\t"+support_after_cutesv_50+"\t"+diff_before_cutesv_50+"\t"+diff_after_cutesv_50+"\n")

    seen_before = get_vcf_distance(args.cutesv_before, truth_large, 500, args.base_bam)
    seen_after = get_vcf_distance(args.cutesv_after, truth_large, 500, args.realigned_bam)
    diff_before_cutesv_500, diff_after_cutesv_500 = get_distance_diff(seen_before, seen_after)
    support_before_cutesv_500, support_after_cutesv_500 = get_support_diff(seen_before, seen_after)
    out_rt.write(support_before_cutesv_500+"\t"+support_after_cutesv_500+"\t"+diff_before_cutesv_500+"\t"+diff_after_cutesv_500+"\n")
