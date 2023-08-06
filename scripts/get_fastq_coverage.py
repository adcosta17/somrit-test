import pysam
import argparse
import sys
import csv
from collections import defaultdict
import gzip
import re
import json
import operator
from intervaltree import Interval, IntervalTree


parser = argparse.ArgumentParser( description='Get metrics for xtea before and after runs')
parser.add_argument('--fastq', required=True)
args = parser.parse_args()

read_lens = defaultdict(int)
count = 0
total = 0
with pysam.FastxFile(args.fastq) as fh:
    for entry in fh:
        count += 1
        total += len(entry.sequence)
        read_lens[len(entry.sequence)] += 1
    # Average Read Length
    print("Average Read Length:")
    print(total/count)
    # Coverage
    print("Coverage:")
    print(total/3000000000)
    # N50
    print("Read N50:")
    sorted_read_lens = dict(sorted(read_lens.items(), key=operator.itemgetter(0),reverse=True))
    count = 0
    n50 = 0
    for key in sorted_read_lens:
        count += (key)*sorted_read_lens[key]
        if count/total > 0.5:
            n50 = key
            break
    print(n50)
