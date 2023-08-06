import pysam
import argparse
import sys
import csv
import re
from collections import defaultdict
import gzip
import random
from intervaltree import Interval, IntervalTree


parser = argparse.ArgumentParser( description='Downsample Fastqs based on fraction')
parser.add_argument('--output-fastq', required=True)
parser.add_argument('--input-fastq', required=True)
parser.add_argument('--input-tsv', required=True)
parser.add_argument('--output-tsv', required=True)
parser.add_argument('--fraction', type=float, required=True)
parser.add_argument('--seed', type=int, required=True)
args = parser.parse_args()

random.seed(args.seed)
# Read in the tsv
tsv_rows = defaultdict(list)
header = ""
with open(args.input_tsv, 'r') as in_tsv:
    count = 0
    for line in in_tsv:
        if count == 0:
            header = line.strip()
            count = 1
            continue
        row = line.strip().split('\t')
        tsv_rows[row[0]].append(row)

with pysam.FastxFile(args.input_fastq) as fin, open(args.output_fastq 'w') as fout, open(args.output_tsv 'w') as out_tsv:
    out_tsv.write(header+"\n")
    for entry in fin:
        if random.uniform(0, 1) < args.fraction:
            fout.write("@"+entry.name+"\n")
            fout.write(entry.sequence+"\n")
            fout.write("+\n")
            fout.write("="*len(entry.sequence)+"\n")
            if entry.name in tsv_rows:
                for row in tsv_rows[entry.name]:
                    out_tsv.write("\t".join(row)+"\n")




