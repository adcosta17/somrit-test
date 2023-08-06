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

tldr_control = defaultdict(IntervalTree)
somrit_control = defaultdict(IntervalTree)

def match(name, truth_reads):
    for read in truth_reads:
        if read in name:
            return read
    return None

def filter_poly_AT(sequence):
    if len(sequence) < 50:
        return "No_Poly_A_Tail:Insert_to_small"
    start = sequence.upper()[:50]
    end = sequence.upper()[len(sequence)-50:]
    pA = "A"*10
    pT = "T"*10
    if pA in start or pT in start or pA in end or pT in end:
        return True
    return False


def same_family(family, annotation):
    if "LINE" in family and "LINE" in annotation:
        return True
    elif "Alu" in family and "Alu" in annotation:
        return True
    elif "SVA" in family and "SVA" in annotation:
        return True
    else:
        return False


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


def get_tldr_test(file, control_samples, parental_calls, insert_positions):
    count = 0
    insert_seen = {}
    seen = defaultdict(IntervalTree)
    with open(file) as in_calls:
        for line in in_calls:
            if count == 0:
                count += 1
                continue
            row = line.strip().split('\t')
            if "PASS" not in line:
                if "LeftFlankSize" in row[24] or "RightFlankSize" in row[24] or "NoFamily" in row[24] or "NonRemappable" in row[24] or "NoTEAlignment" in row[24] or "UnmapCover<0.5" in row[24]:    
                    continue
                if len(row[21]) < 100  or int(row[9]) < 100 or int(row[8])-int(row[7]) < 100 or int(row[15]) < 1 or int(row[16]) > 1 or float(row[11]) < 0.5:
                    continue
            chrom = row[1]
            if chrom == "chr18":
                continue
            start = int(row[2])
            end = int(row[3])
            nearby = insert_positions[chrom][start:end]
            if len(nearby) > 0:
                for item in nearby:
                    insert_seen[chrom+":"+str(item.begin)+"-"+str(item.end)] = 1
            if row[20] == "NA" or not filter_poly_AT(row[21]):
                #TSD
                continue
            # Check samples
            control = False
            for item in row[17].split(','):
                sample = item.split('|')[0]
                if sample in control_samples:
                    control = True
            if control:
                continue
            nearby = insert_positions[chrom][start:end]
            if len(nearby) > 0:
                for item in nearby:
                    insert_seen[chrom+":"+str(item.begin)+"-"+str(item.end)] = 1
            seq_len = int(row[9])
            # Get the samples that are supported by the insertion
            samples = {}
            for item in row[17].split(','):
                samples[item.split('|')[0]] = int(item.split('|')[1])
            parental_nearby = parental_calls[chrom][start:end]
            if len(parental_nearby) > 0:
                continue
            if len(samples) == 1:
                #print(line.strip())
                #print(row[24])
                if int(row[2]) < int(row[3]):
                    seen[row[1]][int(row[2]):int(row[2])+1] = row
                else:
                    seen[row[1]][int(row[3]):int(row[3])+1] = row
    return seen, insert_seen


def get_somrit_test(file, control_samples, parental_calls, insert_positions):
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
                print(count)
            row = line.strip().split('\t')
            if "PASS" not in line or (int(row[13])-int(row[12]))/int(len(row[4])) < 0.5:
                continue
            chrom = row[0]
            if chrom == "chr18":
                continue
            start = int(row[1])
            end = int(row[2])
            if chrom not in seen:
                seen[chrom] = defaultdict(list)
            i = start
            while i < end:
                seen[chrom][i].append(row)
                i += 1
            nearby = insert_positions[chrom][start:end]
            if len(nearby) > 0:
                for item in nearby:
                    insert_seen[chrom+":"+str(item.begin)+"-"+str(item.end)] = 1
    final_counts = defaultdict(IntervalTree)
    count = 0
    with open(file) as in_calls:
        for line in in_calls:
            if count == 0:
                count += 1
                continue
            count += 1
            if count % 100000 == 0:
                print(count)
            row = line.strip().split('\t')
            if "PASS" not in line or "Polymorphic" in line or "No_TSD_Found" in line or "Has_Poly_A_Tail" not in line or (int(row[13])-int(row[12]))/int(len(row[4])) < 0.5:
                continue
            chrom = row[0]
            if chrom == "chr18":
                continue
            start = int(row[1])
            end = int(row[2])
            final_nearby = final_counts[chrom][start-1000:end+1000]
            parental_nearby = parental_calls[chrom][start:end]
            # Check how many samples we see
            nearby = get_nearby(chrom, start-1000, end+1000, seen)
            if len(nearby) == 1:
                control = False
                samples = {}
                reads = row[5].split(',')
                for r in reads:
                    sample = r.split(':')[0]
                    samples[sample] = 1
                    if sample in control_samples:
                        control = True
                if len(samples) > 1 or control:
                    # Have calls that individually are unique to a sample but together are not
                    continue
                if len(final_nearby) > 0 or len(parental_nearby) > 0:
                    continue
                final_counts[chrom][start:end] = [row]
                #print(samples)
                #print(control_samples)
                #print([row])
            else:
                # have multiple insert supporting reads here, generate a list for the position
                ret = []
                samples = {}
                control = False
                for i in nearby:
                    for item in nearby[i]:
                        if "PASS" in item[6]:
                            ret.append(item)
                        reads = item[5].split(',')
                        for r in reads:
                            sample = r.split(':')[0]
                            samples[sample] = 1
                            if sample in control_samples:
                                control = True
                if len(samples) > 1 or control:
                    # Have calls that individually are unique to a sample but together are not
                    continue
                if len(final_nearby) > 0 or len(parental_nearby) > 0:
                    continue
                if len(ret) > 0:
                    final_counts[chrom][start:end] = ret
                #print(ret)
    return final_counts, insert_seen


parser = argparse.ArgumentParser( description='Get metrics for xtea and tldr before and after runs')
parser.add_argument('--tldr', required=True)
parser.add_argument('--somrit', required=True)
parser.add_argument('--insertions', required=True)
parser.add_argument('--control-samples', required=True)
parser.add_argument('--parental', required=True)
args = parser.parse_args()


count = 0
insert_positions_count = 0
insert_positions = defaultdict(IntervalTree)
with open(args.insertions, 'r') as in_csv:
    for line in in_csv:
        if count == 0:
            count = 1
            continue
        row = line.strip().split(',')
        if len(row) > 4:
            if row[2] == "chr18":
                continue
            insert_positions[row[2]][int(row[3])-1000:int(row[4])+1000] = 1
            #print(row[2]+":"+row[3]+"-"+row[4])
            insert_positions_count += 1

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

tldr_calls, tldr_inserts = get_tldr_test(args.tldr, args.control_samples, parental_calls, insert_positions)
somrit_calls, somrit_inserts = get_somrit_test(args.somrit, args.control_samples, parental_calls, insert_positions)

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
        nearby = tldr_calls[chrom][item.begin-1000:item.end+1000]
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

print(insert_positions_count)
print(str(in_set_somrit)+"\t"+str(len(somrit_inserts)))
print(str(in_set_tldr)+"\t"+str(len(tldr_inserts)))
print(somrit_single_read)
