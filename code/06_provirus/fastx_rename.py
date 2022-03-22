#!/usr/bin/env python

import argparse
import os
import re
from Bio import SeqIO
import sys

# create the parser
parser = argparse.ArgumentParser(description='''
   Filter fa file based on sequence lengths
   ''')

# add the arguments
parser.add_argument('--i', help='file to be filtered')
parser.add_argument('--o', help='outfile')
parser.add_argument('--min', help='minimum length', type=int, default=2000)

# parse the command line
args = parser.parse_args()
#args = parser.parse_args('--i graph_prefix_k29.contig --min 100'.split())

# set paths
home = os.getcwd()

# get filenames
file = args.i
path = os.path.split(file)[0]
filename = os.path.split(file)[1]
fileprefix = os.path.splitext(filename)[0]

# parse
outfile = args.o
outhandle = open(outfile , 'w')

pass_count = 0
total_count = 0

for record in SeqIO.parse(open(args.i, 'r'), 'fasta'):
   total_count += 1
   if len(record) >= args.min:
      new_record_name = record.id + '_contig_prophage'
      record.id = new_record_name 
      record.description = ''
      SeqIO.write(record, outhandle, 'fasta')
      pass_count += 1

sys.stderr.write('%s: Written %i of %i (%f)\n' % (outfile, pass_count, total_count, pass_count/total_count*100))

