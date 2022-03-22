#!/bin/python

import sys
import os
import argparse
import csv
import numpy as np
from collections import defaultdict
from collections import Counter

parser = argparse.ArgumentParser(description='''
    Annotating contigs to filter off possible non viral sequences
   ''')
parser.add_argument('-m', help='MAG genome gtdb-tk taxonomy')
parser.add_argument('-o', help='Outfile')

def load_MAG_taxonomy(args):
    '''
    Parse GTDB-TK taxonomy summary file
    '''
    genometax_file = args.m
    lineage = ['superkingdom','phylum','class','order','family','genus','species']
    parsed_motherbins = set()
    genome_tax = dict()
    with open(genometax_file,'r') as csvfile:
        reader = csv.DictReader(csvfile,delimiter='\t')
        for row in reader:
            g = row['user_genome']
            g = g.replace('.fna','')
            tax = row['classification']
            genome_tax[g] = {lin:None for lin in lineage}
            tax = tax.split(';')
            if len(tax) < 2:
                continue
            tax = [ x.split('__')[1] for x in tax ]
            for i,lin in enumerate(lineage):
                entry = tax[i]
                if entry == '':
                    entry = 'NA'
                if '_' in entry and lin != 'species':
                    entry = entry.split('_')[0]
                if '_' in entry and lin == 'species':
                    entries = entry.split(' ')
                    for i,e in enumerate(entries):
                        entries[i] = entries[i].split('_')[0]
                    entry = ' '.join(entries)
                if entry == 'Bacteroidota':
                    entry = 'Bacteroidetes'
                elif entry == 'Verrucomicrobiota':
                    entry = 'Verrucomicrobia'
                elif entry == 'Fusobacteriota':
                    entry = 'Fusobacteria'
                elif entry == 'Synergistetes':
                    entry = 'Synergistetes'
                elif entry == 'Actinobacteriota':
                    entry = 'Actinobacteria'
                elif entry == 'Phocaeicola':
                    entry = 'Bacteroides'
                elif entry == 'Marinifilaceae':
                    entry = 'Odoribacteraceae'
                if lin == 'species' and 'Phocaeicola' in entry:
                    entry = entry.split(' ')
                    entry[0] = 'Bacteroides'
                    entry = ' '.join(entry)
                genome_tax[g][lin] = entry
    return genome_tax


if __name__ == "__main__":
    args = parser.parse_args()

    MAG_tax = load_MAG_taxonomy(args)

    lineage = ['superkingdom','phylum','class','order','family','genus','species']
    with open(args.o,'w') as out:
        header = ['binid']
        for k in lineage:
            header += [k]
        out.write('\t'.join(header)+'\n')
        for binid in MAG_tax.keys():
            if binid == 'user_genome':
                continue
            lineout = [binid]
            for k in lineage:
                if k is None:
                    lineout += ['NA']
                else:
                    lineout += [MAG_tax[binid][k]]
            out.write('\t'.join(lineout)+'\n')

