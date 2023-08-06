import pysam
import argparse
import sys
import csv
import re
from collections import defaultdict
import random
from intervaltree import Interval, IntervalTree

def check_nearby(chrom, start, end, seen):
    i = start
    count = 0
    nearby = {}
    if chrom not in seen:
        return False
    while i < end:
        if i in seen[chrom]:
            return True
        i += 1
    return False

def get_polya_len(seq):
    i = -1
    n = 0 
    while(n < len(seq)):
        if seq[i] != "A":
            break
        i -= 1
        n += 1
    return n

def get_tsd(seq):
    n = random.randint(6,22)
    return seq[len(seq)-n-1:]

def get_insert(seqs, repeat_name):
    n = random.randint(1,2)
    rt = ""
    if n == 1:
        rt = "Alu"
    else:
        rt = "SVA"
    if len(repeat_name) > 0:
        if repeat_name in seqs[rt]:
            sequence = seqs[rt][repeat_name]
            return[repeat_name, sequence, get_polya_len(sequence)]
    n = random.randint(0, len(seqs[rt])-1)
    name = list(seqs[rt].keys())[n]
    sequence = seqs[rt][name]
    return[name, sequence, get_polya_len(sequence)]


def get_insertion_seqs(chromosomes, chromosome_lengths, centromeres, seqs, alu_positions, sva_positions):
    chromosome_list = list(chromosomes.keys())
    n = random.randint(0, len(chromosome_list)-1)
    chrom = chromosome_list[n]
    # Randomly select a position on this chromosome to make an insertion into and select the best contig at the position
    pos = random.randint(50000,chromosome_lengths[chrom]-50001)
    nearby = centromeres[chrom][pos:pos+1]
    nearby_alu = check_nearby(chrom, pos, pos+1, alu_positions)
    nearby_sva = check_nearby(chrom, pos, pos+1, sva_positions)
    found = False
    i = 0
    while(i < 50000):
        j = 0
        while(len(nearby) > 0 and nearby_alu and nearby_sva):
            pos = random.randint(50000,chromosome_lengths[chrom]-50001)
            nearby = centromeres[chrom][pos:pos+1]
            nearby_alu = check_nearby(chrom, pos-10000, pos+10000, alu_positions)
            nearby_sva = check_nearby(chrom, pos-10000, pos+10000, sva_positions)
            j += 1
            if j > 50000:
                break
        insert = get_insert(seqs, "")
        tsd = get_tsd(chromosomes[chrom][pos-50000:pos])
        out_seq = chromosomes[chrom][pos-50000:pos] + insert[1] + tsd + chromosomes[chrom][pos:pos+50000]
        return [out_seq, "BaseRef"+"\tNA\t"+chrom+"\t"+str(pos)+"\t"+insert[0]+"\t"+insert[1]+"\t"+str(insert[2])+"\t"+tsd]
    return ""

def get_insertion_seqs_with_element(chromosomes, chromosome_lengths, centromeres, seqs, alu_positions, sva_positions):
    n = random.randint(1,2)
    rt = ""
    placement = "Start"
    if n == 1:
        rt = "Alu"
    else:
        rt = "SVA"
    chromosome_list = list(chromosomes.keys())
    n = random.randint(0, len(chromosome_list)-1)
    chrom = chromosome_list[n]
    # Randomly select a position on this chromosome to make an insertion into and select the best contig at the position
    repeat_name = ""
    if rt == "Alu":
        positions = list(alu_positions[chrom].keys())
        n = random.randint(0, len(positions)-1)
        pos = positions[n]
        repeat_name = alu_positions[chrom][positions[n]][1]
        m = random.randint(1,2)
        if m == 1:
            pos = alu_positions[chrom][positions[n]][0]
            repeat_name = alu_positions[chrom][positions[n]][1]
            placement = "End"
    else:
        positions = list(sva_positions[chrom].keys())
        n = random.randint(0, len(positions)-1)
        pos = positions[n]
        repeat_name = sva_positions[chrom][positions[n]][1]
        m = random.randint(1,2)
        if m == 1:
            pos = sva_positions[chrom][positions[n]][0]
            placement = "End"
            repeat_name = sva_positions[chrom][positions[n]][1]
    insert = get_insert(seqs, repeat_name)
    tsd = get_tsd(chromosomes[chrom][pos-50000:pos])
    out_seq = chromosomes[chrom][pos-50000:pos] + insert[1] + tsd + chromosomes[chrom][pos:pos+50000]
    return [out_seq, rt+"\t"+placement+"\t"+chrom+"\t"+str(pos)+"\t"+insert[0]+"\t"+insert[1]+"\t"+str(insert[2])+"\t"+tsd]

def new_base(base):
    chars = "ACGT"
    n = random.randint(0,3)
    while base == chars[n]:
        n = random.randint(0,3)
    return chars[n]

def modify_seq(seq, frac):
    # Made random mutations to the sequence at a set probability
    seq_list = list(seq)
    n = int(frac*len(seq_list))
    if n == 0:
        return seq
    seen = {}
    for i in range(n):
        m = random.randint(0, len(seq)-1)
        if m in seen:
            while m in seen:
                m = random.randint(0, len(seq)-1)
        seen[m] = 1
        base = seq_list[m]
        seq_list[m] = new_base(base)
    return "".join(seq_list)

def get_insertion_seqs_duplication(chromosomes, chromosome_lengths, centromeres, seqs, alu_positions, sva_positions, frac):
    n = random.randint(1,2)
    rt = ""
    placement = "Start"
    if n == 1:
        rt = "Alu"
    else:
        rt = "SVA"
    chromosome_list = list(chromosomes.keys())
    n = random.randint(0, len(chromosome_list)-1)
    chrom = chromosome_list[n]
    # Randomly select a position on this chromosome to make an insertion into and select the best contig at the position
    repeat_name = ""
    if rt == "Alu":
        positions = list(alu_positions[chrom].keys())
        n = random.randint(0, len(positions)-1)
        pos_start = positions[n]
        pos_end = pos = alu_positions[chrom][positions[n]][0]
        repeat_name = alu_positions[chrom][positions[n]][1]
        insert = chromosomes[chrom][pos_start:pos_end]
    else:
        positions = list(sva_positions[chrom].keys())
        n = random.randint(0, len(positions)-1)
        pos_start = positions[n]
        pos_end = pos = sva_positions[chrom][positions[n]][0]
        repeat_name = sva_positions[chrom][positions[n]][1]
        insert = chromosomes[chrom][pos_start:pos_end]
    insert = modify_seq(insert, frac)
    out_seq = chromosomes[chrom][pos-50000:pos] + insert + chromosomes[chrom][pos:pos+50000]
    return [out_seq, rt+"Duplication_"+str(frac)+"_\t"+placement+"\t"+chrom+"\t"+str(pos)+"\t"+repeat_name+"\t"+insert+"\tNA\tNA"]


parser = argparse.ArgumentParser( description='Add insertion sequences to contigs')
parser.add_argument('--input', required=True)
parser.add_argument('--output-folder', required=True)
parser.add_argument('--output-prefix', required=True)
parser.add_argument('--output-tsv', required=True)
parser.add_argument('--output-data', required=True)
parser.add_argument('--pbsim-model', required=True)
parser.add_argument('--pbsim-path', required=True)
parser.add_argument('--ref', required=True)
parser.add_argument('--rpmk', required=True)
parser.add_argument('--centromeres', required=True)
parser.add_argument('--telomeres', required=True)
parser.add_argument('--seqs', required=True)
parser.add_argument('--chrom-list', default="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY")
args = parser.parse_args()

telomere_positions= defaultdict(IntervalTree)
centromere_positions = defaultdict(IntervalTree)
with open(args.telomeres) as in_tf:
    count = 0
    for line in in_tf:
        if count == 0:
            count = 1
            continue
        line_args = line.strip().split('\t')
        chrom = line_args[1]
        start = int(line_args[2])
        end = int(line_args[3])
        key = chrom + ":" + str(start) + "-" + str(end)
        telomere_positions[chrom][start:end] = key
with open(args.centromeres) as in_cf:
    count = 0
    for line in in_cf:
        if count == 0:
            count = 1
            continue
        line_args = line.strip().split('\t')
        chrom = line_args[1]
        start = int(line_args[2])
        end = int(line_args[3])
        key = chrom + ":" + str(start) + "-" + str(end)
        centromere_positions[chrom][start:end] = key

seqs = {}
seqs["Alu"] = defaultdict(str)
seqs["SVA"] = defaultdict(str)
with pysam.FastxFile(args.seqs) as fh:
    for entry in fh:
        if "Alu" in entry.name:
            name = entry.name.split('#')[0]
            seqs["Alu"][name] = entry.sequence.upper()
        elif "SVA" in entry.name:
            name = entry.name.split('#')[0]
            seqs["SVA"][name] = entry.sequence.upper()

alu_positions = {}
sva_positions = {}
with open(args.rpmk, 'r') as in_repeat:
    for line in in_repeat:
        row = line.strip().split('\t')
        nearby = centromere_positions[row[5]][int(row[6]):int(row[7])]
        if len(nearby) > 0:
            continue
        if row[12] == "Alu" and row[10] in seqs["Alu"] and (int(row[14])-max(0,int(row[13]))) > 250:            
            if row[5] not in alu_positions:
                alu_positions[row[5]] = {}
            alu_positions[row[5]][int(row[6])] = [int(row[7]), row[10]]
        elif row[12] == "SVA" and row[10] in seqs["SVA"] and (int(row[14])-max(0,int(row[13]))) > 750:            
            if row[5] not in sva_positions:
                sva_positions[row[5]] = {}
            sva_positions[row[5]][int(row[6])] = [int(row[7]), row[10]]


chromosomes = {}
chromosome_lengths = {}
fa = pysam.FastaFile(args.ref)
for chrom in args.chrom_list.split(','):
    seq = fa.fetch(chrom)
    chromosomes[chrom] = seq
    chromosome_lengths[chrom] = len(seq)

count = 0
# Randomly select poisitions from the reference that aren't near centromere or telemere regions
# Check to see if they are within 1000bp of any known SVA or ALU
# If not then insert an SVA or Alu and simulate reads 
with open(args.output_tsv, 'w') as out_tsv, open(args.output_data, 'w') as out_data:
    while count < 1000:
        ret = get_insertion_seqs(chromosomes, chromosome_lengths, centromere_positions, seqs, alu_positions, sva_positions)
        if len(ret) > 0:
            out_data.write(str(count)+"\t"+ret[1]+"\n")
            with open(args.output_folder+"/"+args.output_prefix+"."+str(count)+".fa", 'w') as out_fa:
                out_fa.write(">"+str(count)+"\n"+ret[0]+"\n")
            fa = args.output_folder+"/"+args.output_prefix+"."+str(count)+".fa"
            print(args.pbsim_path+" "+fa+" --prefix "+args.output_folder+"/"+args.output_prefix+"."+str(count)+".pbsim --difference-ratio 23:31:46 --seed 10 --length-min 1000 --depth 100 --hmm_model "+args.pbsim_model+" --length-mean 30000 --accuracy-mean 0.95 --accuracy-min 0.85 --accuracy-max 0.98")
            print("rm "+args.output_folder+"/"+args.output_prefix+"."+str(count)+".pbsim_0001.maf")
            print("rm "+args.output_folder+"/"+args.output_prefix+"."+str(count)+".pbsim_0001.ref")
            out_tsv.write(str(count)+"\t"+args.output_folder+"/"+args.output_prefix+"."+str(count)+".pbsim_0001.fastq\t"+fa+"\n")
            count += 1
    count = 1000
    while count < 2000:
        ret = get_insertion_seqs_with_element(chromosomes, chromosome_lengths, centromere_positions, seqs, alu_positions, sva_positions)
        if len(ret) > 0:
            out_data.write(str(count)+"\t"+ret[1]+"\n")
            with open(args.output_folder+"/"+args.output_prefix+"."+str(count)+"_nested.fa", 'w') as out_fa:
                out_fa.write(">"+str(count)+"\n"+ret[0]+"\n")
            fa = args.output_folder+"/"+args.output_prefix+"."+str(count)+"_nested.fa"
            print(args.pbsim_path+" "+fa+" --prefix "+args.output_folder+"/"+args.output_prefix+"."+str(count)+"_nested.pbsim --difference-ratio 23:31:46 --seed 10 --length-min 1000 --depth 100 --hmm_model "+args.pbsim_model+" --length-mean 30000 --accuracy-mean 0.95 --accuracy-min 0.85 --accuracy-max 0.98")
            print("rm "+args.output_folder+"/"+args.output_prefix+"."+str(count)+"_nested.pbsim_0001.maf")
            print("rm "+args.output_folder+"/"+args.output_prefix+"."+str(count)+"_nested.pbsim_0001.ref")
            out_tsv.write(str(count)+"\t"+args.output_folder+"/"+args.output_prefix+"."+str(count)+"_nested.pbsim_0001.fastq\t"+fa+"\n")
            count += 1
    count = 0
    while count < 5099:
        frac = int(count/100)/100
        ret = get_insertion_seqs_duplication(chromosomes, chromosome_lengths, centromere_positions, seqs, alu_positions, sva_positions, frac)
        if len(ret) > 0:
            out_data.write(str(count+2000)+"\t"+ret[1]+"\n")
            with open(args.output_folder+"/"+args.output_prefix+"."+str(count+2000)+"_duplicate.fa", 'w') as out_fa:
                out_fa.write(">"+str(count+2000)+"\n"+ret[0]+"\n")
            fa = args.output_folder+"/"+args.output_prefix+"."+str(count+2000)+"_duplicate.fa"
            print(args.pbsim_path+" "+fa+" --prefix "+args.output_folder+"/"+args.output_prefix+"."+str(count+2000)+"_duplicate.pbsim --difference-ratio 23:31:46 --seed 10 --length-min 1000 --depth 100 --hmm_model "+args.pbsim_model+" --length-mean 30000 --accuracy-mean 0.95 --accuracy-min 0.85 --accuracy-max 0.98")
            print("rm "+args.output_folder+"/"+args.output_prefix+"."+str(count+2000)+"_duplicate.pbsim_0001.maf")
            print("rm "+args.output_folder+"/"+args.output_prefix+"."+str(count+2000)+"_duplicate.pbsim_0001.ref")
            out_tsv.write(str(count+2000)+"\t"+args.output_folder+"/"+args.output_prefix+"."+str(count+2000)+"_duplicate.pbsim_0001.fastq\t"+fa+"\n")
        count += 1





