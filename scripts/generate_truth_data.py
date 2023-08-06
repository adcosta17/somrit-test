import pysam
import argparse
import sys
import re
from collections import defaultdict

def calculate_sdust_score(seq):
    if seq == "":
        return 0
    triplet_counts = defaultdict(int)
    for i in range(0, len(seq) - 2):
        triplet_counts[seq[i:i+3]] += 1
    sum_score = 0
    for triplet in triplet_counts:
        count = triplet_counts[triplet]
        s = float(count * (count - 1)) / 2.0
        sum_score += s
    if len(seq) - 1 == 0:
        return 0
    sum_score /= (len(seq) - 1)
    return sum_score


def get_deletion_pos(cigarstring):
    # Count up the position on the read until we get to a deletion
    deletion_positions = []
    count = 0
    for cg in re.findall('[0-9]*[A-Z]', cigarstring):
        if cg.endswith('M'):
            count += int(cg[:cg.find("M")])
        elif cg.endswith('X'):
            count += int(cg[:cg.find("X")])
        elif cg.endswith('I'):
            count += int(cg[:cg.find("I")])
        elif cg.endswith('P'):
            count += int(cg[:cg.find("P")])
        elif cg.endswith('S'):
            count += int(cg[:cg.find("S")])
        elif cg.endswith('H'):
            count += int(cg[:cg.find("H")])
        elif cg.endswith('D'):
            if int(cg[:cg.find("D")]) >= 100:
                deletion_positions.append((count,int(cg[:cg.find("D")])))
    return deletion_positions

def update_annotation(annotation, update):
    if annotation == "PASS":
        annotation = update
    else:
        annotation = annotation+","+update
    return annotation


def get_hard_start(cigarstring):
    for cg in re.findall('[0-9]*[A-Z]', cigarstring):
        if cg.endswith('H'):
            return int(cg[:cg.find("H")])
        else:
            break
    return 0

def get_insertion_pos(cigarstring, min_detected_inclusion_length):
    # Count up the position on the read until we get to a deletion
    insert_positions = []
    read_count = 0
    ref_count = 0
    for cg in re.findall('[0-9]*[A-Z=]', cigarstring):
        if cg.endswith('M'):
            read_count += int(cg[:cg.find("M")])
            ref_count += int(cg[:cg.find("M")])
        elif cg.endswith('X'):
            read_count += int(cg[:cg.find("X")])
            ref_count += int(cg[:cg.find("X")])
        elif cg.endswith('='):
            read_count += int(cg[:cg.find("=")])
            ref_count += int(cg[:cg.find("=")])
        elif cg.endswith('I'):
            if int(cg[:cg.find("I")]) > min_detected_inclusion_length:
                insert_positions.append([read_count, int(cg[:cg.find("I")]), ref_count])
            read_count += int(cg[:cg.find("I")])
        elif cg.endswith('S'):
            read_count += int(cg[:cg.find("S")])
        elif cg.endswith('H'):
            read_count += int(cg[:cg.find("H")])
        elif cg.endswith('D'):
            ref_count += int(cg[:cg.find("D")])
        elif cg.endswith('N'):
            ref_count += int(cg[:cg.find("N")])
    return insert_positions

def get_nearby_rm(chrom, start, end, seen):
    ret = {}
    if chrom not in seen:
        return ret
    i = start
    while i <= end:
        if i in seen[chrom]["start"]:
            ret[seen[chrom]["start"][i][1]] = [i, seen[chrom]["start"][i][0]]
        if i in seen[chrom]["end"]:
            ret[seen[chrom]["end"][i][1]] = [seen[chrom]["end"][i][0], i]
        i+=1
    # May be spanning within a block
    for item in seen[chrom]["start"]:
        if item <= start and seen[chrom]["start"][item][0] >= end:
            ret[seen[chrom]["start"][item][1]] = [item, seen[chrom]["start"][item][0]]
        if item >= start and seen[chrom]["start"][item][0] <= end:
            ret[seen[chrom]["start"][item][1]] = [item, seen[chrom]["start"][item][0]]
    return ret

def get_coverage(item, read_start, read_end):
    ann_start = item[0]
    ann_end = item[1]
    if ann_start <= read_start:
        if ann_end <= read_end:
            return ann_end - read_start
        else:
            # Annotation is past the read end
            return read_end - read_start
    else:
        # Annotation starts after the read start
        if ann_end <= read_end:
            return ann_end - ann_start
        else:
            # Annotation is past the read end
            return read_end - ann_start
    return 0

def get_family(item1, item2, item3):
    if "alu" in item1.lower() or "alu" in item2.lower() or "alu" in item3.lower():
        return "alu"
    if "sva" in item1.lower() or "sva" in item2.lower() or "sva" in item3.lower():
        return "sva"
    if "line" in item1.lower() or "line" in item2.lower() or "line" in item3.lower():
        return "line"
    if "simple_repeat" in item1.lower() or "simple_repeat" in item2.lower() or "simple_repeat" in item3.lower():
        return "simple_repeat"
    if "erv" in item1.lower() or "erv" in item2.lower() or "erv" in item3.lower():
        return "erv"
    return item1+"_"+item2+"_"+item3

def get_rm_tagged_inserts(record, rm_list, min_detected_inclusion_length, min_mapq, in_fq):
    records_to_output = {}
    read_annotation = "PASS"
    if record.mapq < min_mapq:
        read_annotation = "mapq<20" 
    # check for any long insertions
    insert_positions = get_insertion_pos(record.cigarstring, min_detected_inclusion_length)
    count = 0
    merged = False
    read_seq = None
    orientation = '+'
    if not record.is_forward:
        orientation = '-'
    if len(insert_positions) > 0 :
        try:
            read_seq = in_fq.fetch(record.query_name)
        except:
            return
    for insert in insert_positions:
        ref_start = insert[2]+record.reference_start
        ref_end = insert[2]+1+record.reference_start
        read_start = insert[0]
        read_end = insert[0]+insert[1]
        if orientation == '-' and read_seq is not None:
            tmp = len(read_seq) - read_start
            read_start = len(read_seq) - read_end
            read_end = tmp
        nearby = get_nearby_rm(record.query_name, read_start, read_end, rm_list)
        rm_tag = ""
        coverage = defaultdict(int)
        for item in nearby:
            cov = get_coverage(nearby[item], read_start, read_end)
            coverage[item] += cov
        for tag in coverage:
            if coverage[tag]/(read_end - read_start) > 0.5:
                if len(rm_tag) > 0:
                    rm_tag += ","
                rm_tag += tag
        if orientation == '-' and read_seq is not None:
            tmp = len(read_seq) - read_start
            read_start = len(read_seq) - read_end
            read_end = tmp
        insertion_sequence = ""
        if read_seq is not None:
            insertion_sequence = read_seq[read_start:read_end]
        annotation = read_annotation
        output = "%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%s\t%s" % (record.reference_name, ref_start, ref_end, record.query_name, read_start, read_end, orientation, insertion_sequence,annotation, rm_tag)
        print(output)
        count += 1

parser = argparse.ArgumentParser( description='Extract reads with long insertions')
parser.add_argument('--pat-bam', type=str, required=True)
parser.add_argument('--mat-bam', type=str, required=True)
parser.add_argument('--pat-fasta', type=str, required=True)
parser.add_argument('--mat-fasta', type=str, required=True)
parser.add_argument('--min-insertion-length', type=int, default=100)
parser.add_argument('--min-detected-inclusion-length', type=int, default=100)
parser.add_argument('--min-flank-size', required=False, default=100)
parser.add_argument('--min-mapq', required=False, type=int, default=20)
parser.add_argument('--reference-gap-minimum', type=int, default=100)
parser.add_argument('--repeat-master-pat', type=str, required=True)
parser.add_argument('--repeat-master-mat', type=str, required=True)
args = parser.parse_args()

print("\t".join(["chromosome", "reference_insertion_start", "reference_insertion_end", "read_name", "read_insertion_start", "read_insertion_end", "dust_score", "insertion_sequence","pass_fail", "annotation"]))

# mat
mat_rm = {}
with open(args.repeat_master_mat, 'r') as in_rm:
    for line in in_rm:
        row = line.strip().split('\t')
        if row[0] not in mat_rm:
            mat_rm[row[0]] = {}
            mat_rm[row[0]]["start"] = {}
            mat_rm[row[0]]["end"] = {}
        mat_rm[row[0]]["start"][int(row[1])] = [int(row[2]), get_family(row[3], row[6], row[7])]
        mat_rm[row[0]]["end"][int(row[2])] = [int(row[1]), get_family(row[3], row[6], row[7])]
in_fq = pysam.FastaFile(args.mat_fasta)
sam_reader = pysam.AlignmentFile(args.mat_bam)
for record in sam_reader.fetch():
    get_rm_tagged_inserts(record, mat_rm, args.min_detected_inclusion_length, args.min_mapq, in_fq)


#pat
pat_rm = {}
with open(args.repeat_master_pat, 'r') as in_rm:
    for line in in_rm:
        row = line.strip().split('\t')
        if row[0] not in pat_rm:
            pat_rm[row[0]] = defaultdict(str)
            pat_rm[row[0]]["start"] = {}
            pat_rm[row[0]]["end"] = {}
        pat_rm[row[0]]["start"][int(row[1])] = [int(row[2]), get_family(row[3], row[6], row[7])]
        pat_rm[row[0]]["end"][int(row[2])] = [int(row[1]), get_family(row[3], row[6], row[7])]
in_fq = pysam.FastaFile(args.pat_fasta)
sam_reader = pysam.AlignmentFile(args.pat_bam)
for record in sam_reader.fetch():
    get_rm_tagged_inserts(record, pat_rm, args.min_detected_inclusion_length, args.min_mapq, in_fq)





