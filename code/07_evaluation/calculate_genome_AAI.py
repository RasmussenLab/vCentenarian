#!/usr/bin/env python

import os, numpy as np, Bio.SeqIO, argparse, time

parser = argparse.ArgumentParser()
parser.add_argument('--in_faa', required=True, metavar='PATH')
parser.add_argument('--in_blast', required=True, metavar='PATH')
parser.add_argument('--out_tsv', required=True, metavar='PATH')
args = vars(parser.parse_args())


def get_genome_AAI(outhandle,query,genome_size,genome_alns):
    '''Calculates pairwise AAI for all protein hits'''
    query_genes = genome_size[query]
    for target in genome_alns[query]:
        target_genes = genome_size[target]
        alns = genome_alns[query][target].values()
        shared_genes = len(alns)
        qcov = 100.0*shared_genes/query_genes
        tcov = 100.0*shared_genes/target_genes
        aai = np.mean([float(_[2]) for _ in alns])
        row = [query, target, query_genes, target_genes, shared_genes, qcov, tcov, aai]
        outhandle.write('\t'.join([str(_) for _ in row])+'\n')


from collections import defaultdict
genome_size = defaultdict(int)
gene_to_genome = {}
for index, r in enumerate(Bio.SeqIO.parse(args['in_faa'], 'fasta')):
	genome_id = r.id.rsplit('_', 1)[0]
	genome_size[genome_id] += 1
	gene_to_genome[r.id] = genome_id

genome_alns = {}
current_genome = None
print('Parsing blastfile...')
start = time.time()

out = open(args['out_tsv'], 'w')
header = ['qname', 'tname', 'qgenes', 'tgenes', 'sgenes', 'qcov', 'tcov', 'aai']
out.write('\t'.join(header)+'\n')
for index, line in enumerate(open(args['in_blast'])):
    aln = line.split()
    gene = aln[0]
    query = gene_to_genome[aln[0]]
    target = gene_to_genome[aln[1]]
    score = float(aln[-1])

    if query != current_genome:
        if not current_genome is None:
            get_genome_AAI(out, current_genome,genome_size, genome_alns)
            genome_alns[query] = {}
            current_genome = query
        else:
            genome_alns[query] = {}
            current_genome = query
    if target not in genome_alns[query]:
        genome_alns[query][target] = {}
    if gene not in genome_alns[query][target]:
        genome_alns[query][target][gene] = aln
    elif score > float(genome_alns[query][target][gene][-1]):
        genome_alns[query][target][gene] = aln

### The last genome processed
get_genome_AAI(out, current_genome,genome_size, genome_alns)
out.close()

end = time.time()
print(end - start)
