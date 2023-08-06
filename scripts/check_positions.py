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
parser.add_argument('--inserts-tsv', required=True)
parser.add_argument('--bam', required=True)
args = parser.parse_args()

# Start by reading in the list of spiked in reads and the positions they come from
# Then read in the metadata on those positions that will tell us more about the positions
insert_data = defaultdict(list)

with open(args.inserts_tsv, 'r') as in_tsv:
    for line in in_tsv:
        row = line.strip().split('\t')
        insert_data[row[0]] = row

samfile = pysam.AlignmentFile(args.bam)
seen = defaultdict(list)
for record in samfile.fetch():
    if record.reference_name == insert_data[record.query_name][5] and record.reference_start < int(insert_data[record.query_name][6]) and record.reference_end > int(insert_data[record.query_name][6]):
        # have a matching hit
        seen[record.query_name].append("Match\t"+record.reference_name+"\t"+str(record.reference_start)+"\t"+str(record.reference_end)+"\t"+str(record.mapping_quality)+"\t"+"\t".join(insert_data[record.query_name]))
    else:
        seen[record.query_name].append("Elsewhere\t"+record.reference_name+"\t"+str(record.reference_start)+"\t"+str(record.reference_end)+"\t"+str(record.mapping_quality)+"\t"+"\t".join(insert_data[record.query_name]))

for pos in insert_data:
    if pos not in seen:
        seen[record.query_name].append("Missing\t"+record.reference_name+"\t"+str(record.reference_start)+"\t"+str(record.reference_end)+"\t"+str(record.mapping_quality)+"\t"+"\t".join(insert_data[record.query_name]))

for pos in seen:
    for item in seen[pos]:
        print(item)

