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
parser.add_argument('-g', help='Genome meta-taxonomy file')
parser.add_argument('-m', help='MAG genome gtdb-tk taxonomy')
parser.add_argument('-b', help='Spacer blast file for MAGs')
parser.add_argument('-c', help='Spacer blast file for NCBI references')
parser.add_argument('-v', help='Viruses in db')
parser.add_argument('-o', help='Outfile')


class Virus_class:
    def __init__(self):
        self.domintant_spacer_lineage = None
        self.domintant_prophage_lineage = None
        self.purity_spacer_lineage = None
        self.purity_prophage_lineage = None

def load_genome_taxonomy(args):
    '''
    Loads Genome taxonomy associated with Each Spacer in the Spacer-database
    '''
    spacertax_file = args.g

    genome_tax = dict()
    with open(spacertax_file,'r') as csvfile:
        reader = csv.DictReader(csvfile,delimiter='\t')
        for row in reader:
            organism = row['spacerid']
            genome_tax[organism] = row

    return genome_tax

def load_MAG_taxonomy(args):
    '''
    Parse GTDB-TK taxonomy summary file
    '''
    genometax_file = args.m
    lineage = ['superkingdom','phylum','class','order','family','genus','species']
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


def parse_CRISPR_results(args):
    '''
    NOT used anymore
    The CRISPR-open database is rather redundant. Do not count a spacer-hit from the same genome more than once. 
    '''
    ### Parse blastn 
    virus_genome_dict = {}
    virus_host = dict()
    spaceromefile = args.c
    with open(spaceromefile,'r') as infile:
        for line in infile:
            line = line.strip().split()
            spacername, virus, identity, spacer_covered, spacer_len = line[0], line[1], float(line[2]), int(line[3]),int(line[-2])
            genome = spacername.split('.')[0]
            bitscore = float(line[-3])
            evalue = float(line[-4])
            mismatches = int(line[4])
            
            ### Stringent filtering of hits- Do not count spacers from the same genome more than once.
            if mismatches > 2:
                continue
            if identity >= 80 and spacer_covered/spacer_len >=0.90:
                if not virus in virus_host:
                    virus_host[virus] = [spacername]
                    virus_genome_dict[virus] = dict()
                    virus_genome_dict[virus][genome] = ''
                else:
                    if not genome in virus_genome_dict[virus]:
                        virus_host[virus] += [spacername]
    return virus_host

def parse_CRISPR_results_MAG(args):
    '''Parse blast hits towards MAG
        1. Maximum mismatches 2 
        2. Accepted hits spacer coverage > 90% and >80% nucleotide identity 
    '''
    ### Parse blastn 
    virus_host = dict()
    spaceromefile = args.b
    with open(spaceromefile,'r') as infile:
        for line in infile:
            line = line.strip().split()
            spacername, virus, identity, spacer_covered, spacer_len = line[0], line[1], float(line[2]), int(line[3]),int(line[-2])
            genome = '::'.join(spacername.split('::')[:2])
            bitscore = float(line[-3])
            evalue = float(line[-4])
            mismatches = int(line[4])
            
            ### Stringent filtering of hits
            if mismatches > 2:
                continue
            if identity >= 80 and spacer_covered/spacer_len >=0.90:
                if not virus in virus_host:
                    virus_host[virus] = defaultdict(list)
                    virus_host[virus][genome] = [spacername]
                else:
                    virus_host[virus][genome] += [spacername]
    return virus_host


def evaluate_host_consistency(args,viruses):
    """
    2. For each Virus load all Hosts by Spacers 
    4. For both Methods calculate the most common host at each Lineage 
    5. For both Methods calculate Host-assignment overall Purity
    """ 
    lineage = ['superkingdom','phylum','class','order','family','genus','species']
    ### Calculate the most common Viral host across lineages

    class Virus_class:
        def __init__(self):
            self.domintant_spacer_lineage = None
            self.domintant_prophage_lineage = None
            self.purity_spacer_lineage = None
            self.purity_prophage_lineage = None
    

    MAG_tax = load_MAG_taxonomy(args)
    #genome_tax = load_genome_taxonomy(args)

    #virus_host = parse_CRISPR_results(args)
    virus_MAG = parse_CRISPR_results_MAG(args)

    virus_annotations = dict()
    for virus in viruses:
        if not virus in virus_annotations:
            viral_class = Virus_class()
            viral_class.spacers = {k:[] for k in lineage}
        else:
            viral_class = virus_annotations[virus] 

        ### Reference database spacers
        #if virus in virus_host:
        #    for spacer in virus_host[virus]:
        #        organism = spacer.split('.')[0]
        #        glineage = genome_tax[organism]
        #        for k in lineage:
        #            lin = glineage[k].replace('[','').replace(']','')
        #            viral_class.spacers[k].append(lin) 
        
        ### MAG database spacers
        if virus in virus_MAG:
            for g in virus_MAG[virus]:
                
                ### In case the MAG is not annotated
                if not g in MAG_tax:
                    continue
                glineage = MAG_tax[g]
                for i in virus_MAG[virus][g]:
                    for k in lineage:
                        lin = glineage[k].replace('[','').replace(']','')
                        viral_class.spacers[k].append(lin) 
        virus_annotations[virus] = viral_class
            
    ### Calculate most Common Host Taxonomy at each Lineage 
    ### & Evaluate Host-prediction-purity 
    for virus in virus_annotations:
        viral_class = virus_annotations[virus]
        viral_class.domintant_spacer_lineage = None
        viral_class.domintant_prophage_lineage = None
        dominant_spacer_lineage = {k:'NA' for k in lineage}
        purity_spacer_lineage = {k:'NA' for k in lineage}

        for k in lineage:
            c = Counter(viral_class.spacers[k])            
            for entry in c:
                proportion = c[entry]/len(viral_class.spacers[k])
            cmc = c.most_common()
            if len(cmc) != 0:
                purity = viral_class.spacers[k].count(cmc[0][0])/len(viral_class.spacers[k])
                if purity >= 0.4:
                    dominant_spacer_lineage[k] = cmc[0][0]
                purity_spacer_lineage[k] = purity

        if dominant_spacer_lineage['superkingdom'] == 'NA':
            viral_class.domintant_spacer_lineage = None
        else:
            viral_class.domintant_spacer_lineage = dominant_spacer_lineage
            viral_class.purity_spacer_lineage = purity_spacer_lineage

        virus_annotations[virus] = viral_class
    return virus_annotations

if __name__ == "__main__":
    args = parser.parse_args()

    viruses = list()
    with open(args.v, 'r') as infile:
        for line in infile:
            viruses += [line.strip().replace('>','')]
    virus_host_annotations = evaluate_host_consistency(args,viruses)

    ### Write consus results for each virus
    lineage = ['superkingdom','phylum','class','order','family','genus','species']
    with open(args.o,'w') as out:
        header = ['virus']
        for k in lineage:
            header += [k]
        out.write('\t'.join(header)+'\n')
        for virus in viruses:
            lineout = [virus]

            ### No spacers 
            if not virus in virus_host_annotations:
                for k in lineage:
                    lineout += ['NA']
                out.write('\t'.join(lineout)+'\n')
            else:
                ### Some spacers
                viral_class = virus_host_annotations[virus]
                if viral_class.domintant_spacer_lineage is None:
                    for k in lineage:
                        lineout += ['NA']
                    out.write('\t'.join(lineout)+'\n')
                else:
                    dominant_spacer_lineage = viral_class.domintant_spacer_lineage
                    purity_spacer_lineage = viral_class.purity_spacer_lineage
                    for k in lineage:
                        if k == 'species':
                            species = dominant_spacer_lineage[k]
                            genera = species.split(' ')[0]
                            lineout[-1] = genera
                            lineout += [species]
                        else:
                            lineout += [dominant_spacer_lineage[k]]
                    out.write('\t'.join(lineout)+'\n')