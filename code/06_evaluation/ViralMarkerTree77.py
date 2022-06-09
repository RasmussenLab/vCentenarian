#!/bin/python

import os
import Bio.SeqIO
import subprocess as sp
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--in_faa', required=True, metavar='PATH')
parser.add_argument('--out_dir', required=True, metavar='PATH')
parser.add_argument('--threads', type=int, default=1, metavar='INT')
args = vars(parser.parse_args())

# make dir
if not os.path.exists(args['out_dir']):
	os.makedirs(args['out_dir'])

# read proteins
faa = dict([(r.id, str(r.seq).rstrip('*')) for r in Bio.SeqIO.parse(args['in_faa'], 'fasta')])

print ("search 77 markers versus original proteins")
if not os.path.exists(os.path.join(args['out_dir'],'hmmsearch.txt')):
	cmd = "hmmsearch --noali --cpu %s " % args['threads']
	cmd += "--tblout %s/hmmsearch.txt " % args['out_dir']
	cmd += "VOG77.hmm "
	cmd += "%s > /dev/null" % args['in_faa']
	p = sp.Popen(cmd, shell=True)
	p.wait()

print ("extract proteins for top hmm hits")
hits_dir = '%s/hits' % args['out_dir']

if not os.path.exists(hits_dir):
	os.makedirs(hits_dir)
hits = {}
for line in open('%s/hmmsearch.txt' % args['out_dir']):
	if line[0]=='#': continue
	r = line.split()
	if float(r[4]) > 1e-5: continue
	elif r[0] not in hits:
		hits[r[0]] = r
	elif float(r[5]) > float(hits[r[0]][5]):
		hits[r[0]] = r

marker_to_hits = {}
for r in hits.values():
	marker_id = r[2]
	if marker_id not in marker_to_hits:
		marker_to_hits[marker_id] = []
	gene = faa[r[0]]
	marker_to_hits[marker_id].append([r[0], gene])
for marker_id in marker_to_hits:
	outpath = '%s/%s.faa' % (hits_dir, marker_id)
	out = open(outpath, 'w')
	for gene_id, seq in marker_to_hits[marker_id]:
		out.write('>'+gene_id+'\n'+seq+'\n')

print ("build multiple seq alignment (MSA) for each PC, round 2")
msa2_dir = '%s/msa2' % args['out_dir']
if not os.path.exists(msa2_dir):
	os.makedirs(msa2_dir)
for file in os.listdir(hits_dir):
	inpath = hits_dir+'/'+file
	outpath = msa2_dir+'/'+file
	cmd = "famsa -t %s %s %s" % (args['threads'], inpath, outpath)
	p = sp.Popen(cmd, shell=True)
	p.wait()

print ("trim multiple seq alignments")
trim2_dir = '%s/trim2' % args['out_dir']
if not os.path.exists(trim2_dir):
	os.makedirs(trim2_dir)
for file in os.listdir(msa2_dir):
	inpath = msa2_dir+'/'+file
	outpath = trim2_dir+'/'+file
	cmd = "trimal -keepheader -gt 0.5 -in %s -out %s" % (inpath, outpath)
	p = sp.Popen(cmd, shell=True)
	p.wait()

print ("concatenate alignments and fill missing genes with gaps")
genome_ids = set([_[1:].split()[0].rsplit('_', 1)[0] for _ in open(args['in_faa']) if _[0]=='>'])

msa_catted = dict([(id, '') for id in genome_ids])
for index, file in enumerate(os.listdir(trim2_dir)):
	marker_id = file.split('.')[0]
	genome_to_gene = {}
	inpath = '%s/%s' % (trim2_dir, file)
	for r in Bio.SeqIO.parse(inpath, format='fasta'):
		genome_id = r.id.rsplit('_', 1)[0]
		genome_to_gene[genome_id] = r
	marker_length = len(list(genome_to_gene.values())[0].seq)
	for genome_id in genome_ids:
		if genome_id not in genome_to_gene:
			msa_catted[genome_id] += "-" * marker_length
		else:
			gene = genome_to_gene[genome_id]
			msa_catted[genome_id] += str(gene.seq).replace('X', '-')

print ("count gaps per genome")
genome_to_gaps = {}
for genome_id in genome_ids:
	seq = msa_catted[genome_id]
	perc_gaps = round(100*seq.count('-')/float(len(seq)), 2)
	genome_to_gaps[genome_id] = perc_gaps

print ("write concatenated alignment for genomes with <95% gaps ")
with open('%s/concat.faa' % args['out_dir'], 'w') as f:
	for genome_id, seq in msa_catted.items():
		perc_gaps = genome_to_gaps[genome_id]
		if perc_gaps < 95:
			f.write('>%s percent_gaps=%s\n' % (genome_id, perc_gaps))
			f.write('%s\n' % seq)


'''
Then run a tool like IQtree or similar to make the phylogenetic Tree
The CCP77 reference tree was inferred from the MSA using IQ-TREE version 1.5.5 (ref. 51) under the partitioned52 
LG + GAMMA model and midpoint rooted, followed by 100 non-parametric bootstrap replicates to estimate branch support. 

'''