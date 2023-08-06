import pysam
import argparse
import sys
import csv
import re
from collections import defaultdict
import random
from intervaltree import Interval, IntervalTree
import mappy as mp
from multiprocessing import Lock, Queue
import threading
import time
import queue

def get_insertion_pos(cigarstring, ref_start):
    # Count up the position on the read until we get to a deletion
    insert_positions = []
    read_count = 0
    ref_count = ref_start
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
            if int(cg[:cg.find("I")]) > 50:
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

def get_name(n):
    tmp = hex(n)[2:].zfill(32)
    return tmp[0:8]+"-"+tmp[8:12]+"-"+tmp[12:16]+"-"+tmp[16:20]+"-"+tmp[20:]


class mappingThread(threading.Thread):
    def __init__(self, i, q, ref_aligner, insert_dict, t_lock, out_tsv):
        threading.Thread.__init__(self)
        self.i = i
        self.q = q
        self.ref_aligner = ref_aligner
        self.insert_dict = insert_dict
        self.out_str = ""
        self.t_lock = t_lock
        self.out_tsv = out_tsv
    def run(self):
        while(True):
            if self.q.empty():
                print("Done "+str(self.i))
                self.t_lock.acquire()
                self.out_tsv.write(self.out_str)
                self.t_lock.release()
                return
            try:
                data = self.q.get(True, 100)
            except queue.Empty:
                print("Done "+str(self.i))
                self.t_lock.acquire()
                self.out_tsv.write(self.out_str)
                self.t_lock.release()
                return
            row = data
            print(row[0]+"\t"+str(self.i))
            chrom = self.insert_dict[row[0]][3]
            start = int(self.insert_dict[row[0]][4])
            # Map the reads for this sequence to the ref
            nearest = defaultdict(int)
            for name, seq, qual in mp.fastx_read(row[1]): # read a fasta/q sequence
                for hit in self.ref_aligner.map(seq): # traverse alignments
                    if hit.ctg == chrom and hit.mapq > 20 and hit.r_st < (start - 1000) and hit.r_en > (start + 1000):
                        #print(chrom + ":"+ str(start))
                        #print(hit)
                        # Read spans insert pos
                        insert_positions = get_insertion_pos(hit.cigar_str, hit.r_st)
                        for insert in insert_positions:
                            pos = insert[2]
                            #print("Insert at: "+str(pos))
                            if name not in nearest:
                                nearest[name] = abs(pos - start)
                            if abs(pos - start) < nearest[name]:
                                nearest[name] = abs(pos - start)
            within_1000 = 0
            within_500 = 0
            within_250 = 0
            within_100 = 0
            within_50 = 0
            for read in nearest:
                if nearest[read] < 1000:
                    within_1000 += 1
                if nearest[read] < 500:
                    within_500 += 1
                if nearest[read] < 250:
                    within_250 += 1
                if nearest[read] < 100:
                    within_100 += 1
                if nearest[read] < 50:
                    within_50 += 1
            self.out_str += row[0]+"\t"+self.insert_dict[row[0]][1]+"\t"+self.insert_dict[row[0]][2]+"\t"+self.insert_dict[row[0]][5]+"\t"+str(len(nearest))+"\t"+str(within_1000)+"\t"+str(within_500)+"\t"+str(within_250)+"\t"+str(within_100)+"\t"+str(within_50)+"\n"



parser = argparse.ArgumentParser( description='Get reads that span the insertion sequence')
parser.add_argument('--input', required=True)
parser.add_argument('--inserts', required=True)
parser.add_argument('--tsv', required=True)
parser.add_argument('--ref', required=True)
parser.add_argument('--threads', required=False, type=int, default=1)
args = parser.parse_args()

insert_dict = {}
with open(args.inserts, 'r') as in_tsv:
    for line in in_tsv:
        row = line.strip().split('\t')
        insert_dict[row[0]] = row

q = queue.Queue()
with open(args.input, 'r') as in_tsv:
    count = 1
    for line in in_tsv:
        row = line.strip().split('\t')
        q.put(row)

ref_aligner = mp.Aligner(args.ref)
t_lock = threading.Lock()
thread_list = [None] *args.threads
with open(args.tsv, 'w') as out_tsv:
    for i in range(args.threads):
        thread_list[i] = mappingThread( i, q, ref_aligner, insert_dict, t_lock, out_tsv)
        thread_list[i].start()
    print("Setup Threads")
    for i in range(args.threads):
        thread_list[i].join()
