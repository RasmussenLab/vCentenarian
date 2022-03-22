#!/bin/python
import os
from os import listdir
import sys
import csv
import gzip
import numpy as np
from optparse import OptionParser


def load_genome_taxonomy(options):
    '''
    Parse GTDB-TK taxonomy summary file
    '''
    genometax_file = options.genometax
    lineage = ['superkingdom','phylum','class','order','family','genus','species']
    #dlineage = {lin:None for lin in lineage}

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

def load_contig2bin(representative_bins_dir):
    '''Organise contigs in bins'''
    bins = [f for f in listdir(representative_bins_dir)] 
    bin_files = [ os.path.join(representative_bins_dir,b) for b in bins]
    contig2bin = dict() 
    for i, b in enumerate(bin_files):
        binid = bins[i].split('.')[0]
        with gzip.open(b,'rt') as infile:
            for line in infile:
                if line[0] == '>':
                    contig = line.strip().split(' ')[0].replace('>','')
                    contig2bin[contig] = binid
    return contig2bin

def MSP_hits(msp,mspdir,contig2bin):
    '''Parse MSP core-gene hits to MAG bins'''
    blastn_file = os.path.join(mspdir,msp+".core_genes.MAG.m6")
    if os.path.exists(blastn_file):
        binstats = dict()
        core_genes = set()
        gene2len = dict()
        with open(blastn_file,'r') as infile:
            for line in infile:
                line = line.strip().split()
                core_gene , contig, seq_identity, gene_aln, bitscore, core_gene_len, contig_len = line[0], line[1], line[2], line[3], line[-3], line[-2], line[-1]
                core_genes.add(core_gene)
                gene2len[core_gene] = int(core_gene_len)                
                binid = contig2bin[contig]
                if not binid in binstats:
                    binstats[binid] = {'identityxaln':0, 'bitscore':0,'alignment':0,'genes':[]}

                binstats[binid]['identityxaln'] += float(seq_identity)*int(gene_aln)
                binstats[binid]['bitscore'] += float(bitscore)
                binstats[binid]['alignment'] += int(gene_aln)
                binstats[binid]['genes'] += [core_gene]

        ncore_genes = len(core_genes)
        total_core_len = np.sum([gene2len[c] for c in core_genes])
        binsummary = dict() 
        for k,v in binstats.items():
            AF = v['alignment']/total_core_len
            ANI = v['identityxaln']/v['alignment']
            frac_genes = len(set(v['genes']))/ncore_genes
            binsummary[k] = {'AF':AF,'ANI':ANI,'frac_genes':frac_genes,'total_bitscore':v['bitscore']}

        ### Select top3 
        sorted_binsummary = {k: v for k, v in sorted(binsummary.items(), key=lambda value: value[1]['total_bitscore'],reverse=True)}
        top3_bins = list(sorted_binsummary)[:3]
        top3 = {msp:dict()}
        for k in top3_bins:
            top3[msp][k] = binsummary[k]
        return top3


def main():

  parser = OptionParser()
  parser.add_option("-i", "--bindirectory", dest="bindir", help="directory of MAG bins", metavar="PATH")
  parser.add_option("-m", "--msps", dest="mspdir", help="MSPS directory", metavar="PATH")
  parser.add_option("-o", "--outfile", dest="out_file", help="out file", metavar="FILENAME")
  parser.add_option('-g', "--genometax",dest="genometax",help='GTDBTK-taxonomy file')
  (options, args) = parser.parse_args()

  bindir = options.bindir
  mspdir = options.mspdir
  outfile = options.out_file

  ### Pars Bin-fasta files to make a lookup dictionary
  contig2bin = load_contig2bin(bindir)
  genome_tax = load_genome_taxonomy(options)
  lineage = ['superkingdom','phylum','class','order','family','genus','species']
  ### 
  msps = set([f.split('.')[0] for f in listdir(mspdir) if 'msp' in f])
  with open(outfile,'w') as out:
    header = ['msp','binid','AF','ANI','frac_genes','bitscore']
    header = header + lineage
    out.write('\t'.join(header)+'\n')
    for msp in sorted(msps):
        hits =  MSP_hits(msp,mspdir,contig2bin)
        if hits is None:
            NA_line = ['NA' for i in range(len(header)) -1 ]
            lineout = [msp] + NA_line
            out.write('\t'.join(lineout)+'\n')
        else:
            for k,v in hits[msp].items():
                tax = []
                for lin in lineage:
                    tax += [genome_tax[k][lin]]
                lineout = [msp,k, round(v['AF'],2), round(v['ANI'],2), round(v['frac_genes'],2), round(v['total_bitscore'],2)]
                lineout = lineout + tax
                lineout = [str(i) for i in lineout]
                out.write('\t'.join(lineout)+'\n')

if __name__ == '__main__':
  main()



