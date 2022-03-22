#!/bin/python

from collections import defaultdict
from collections import OrderedDict
import collections
from doctest import master
import gzip
import os
import csv
import tree_functions


'''
1. Master-table (Prepared in R...) will contain Genome quality information and an indicator if a vOTU is virus based on 
    1.1. If CheckV quality is Complete (DTR-based), HQ or MQ 
    1.2. If the MGV evaluation pipeline determines the vOTU as virus
    1.3. If >40% of genes are viral annotated or <10% host annotated

2. Add Species-cluster information based on Whole-genome clustering with MGV genomes (ANI-based) - Transfer taxonomy if clustering on species level 
3. Add Genus-cluster information based on Whole-proteoe clustering with MGV genomes (AAI-based) - Transfer genus taxonomy if clustering genus level 

'''

class ViralCluster:
    def __init__(self):
        self.representative = None
        self.taxonomy = None
        self.nmembers = None
        self.taxonomy = None
        self.lifestyle = None
        self.hosttax = None
        self.MGVcluster = False

def VC_taxonomy(rep,mgv_members,master_table):
    tax = None
    if 'MGV' in rep:
        entry = master_table[rep]
        tax = dict((k, entry[k]) for k in ('baltimore','ictv_order','ictv_family','ictv_genus'))
    else:
        entry = master_table[mgv_members[0]]
        tax = dict((k, entry[k]) for k in ('baltimore','ictv_order','ictv_family','ictv_genus'))
    return tax
        
def VC_host(rep,mgv_members,master_table):
    host = None
    if 'MGV' in rep and rep in master_table:
        host = master_table[rep]['host']
    else:
        for cluster_member in mgv_members:
            if cluster_member in master_table and 'MGV' in cluster_member:
                host = master_table[cluster_member]['host']
    return host

def load_species_clusters(master_table, species_file):
    '''Load virus-clustering file'''
    species_file = 'virus.clusters.95'
    virus_VC_lookup = dict()
    VC = dict()
    number_vc = 0
    with open(species_file,'r') as infile:
        for line in infile:
            number_vc += 1
            vc = ViralCluster()
            cluster_representative, members = line.strip().split()
            vc_members = members.split(',')
            mgv_members = []
            for i, cluster_member in enumerate(vc_members):
                if 'MGV' in cluster_member:
                    vc.MGVcluster = True
                    mgv_members += [cluster_member]
                else:
                    virus_VC_lookup[cluster_member] = number_vc
            vc.representative = cluster_representative 
            vc.nmembers = len(vc_members)
            if len(mgv_members) != 0: 
                ### Assign Taxonomy based on MGV's in the CLuster
                vc.taxonomy = VC_taxonomy(cluster_representative, mgv_members, master_table)
                vc.hosttax = VC_host(cluster_representative,mgv_members,master_table)
            VC[number_vc] = vc

    ### Assign VC-information to each vOTU
    tax_entries = ['baltimore','ictv_order','ictv_family','ictv_genus']
    for vOTU in master_table:
        if not 'MGV' in vOTU:
            VC_number = virus_VC_lookup[vOTU]
            VC_name = 'VC'+str(VC_number)
            vc = VC[VC_number]
            if not vc.taxonomy is None: 
                for t in tax_entries:
                    master_table[vOTU][t] = vc.taxonomy[t]
            else:
                for t in tax_entries:
                    master_table[vOTU][t] = 'NA'
            if vc.MGVcluster is True: master_table[vOTU]['MGVref'] = 'MGVref'
            else:
                master_table[vOTU]['MGVref'] = 'Novel'
            master_table[vOTU]['VC'] = VC_name
    return VC, virus_VC_lookup

def load_genus_clusters(master_table,genus_file):
    '''Parse and annotate genus vOTUs in master table'''
    viral_genus_clusters = dict()
    genus_count = 0
    with open(genus_file,'r') as infile:
        for line in infile:
            genus_count += 1
            genus_id = 'G'+str(genus_count)
            vOTUs = line.strip().split()

            if not genus_id in viral_genus_clusters:
                viral_genus_clusters[genus_id] = []

            for virus in vOTUs:
                viral_genus_clusters[genus_id].append(virus)
                if virus in master_table:
                    master_table[virus]['genus'] = genus_id
                else:
                    master_table[virus]['genus'] = 'NA'
    return viral_genus_clusters

def write_VOG_table(master_table,VOGs,votu_markers):
    '''Write table with VOG counts for each vOTU'''
    votu_information = ['contig_id','ictv_order','ictv_family','ictv_genus','genus','decision','length','MGVref']
    header = ['contig_id','ictv_order','ictv_family','ictv_genus','genus','decision','length','MGVref']
    for v in sorted(VOGs): header += [v] 
    header += ['nVOG']   
    with open('VOGtable.tsv','w') as out:
        out.write('\t'.join(header)+'\n')
        for votu in master_table:
            row = []
            i = 0 
            for k in votu_information: row += [master_table[votu][k]]
            for v in sorted(VOGs):
                if votu in votu_markers:
                    if v in votu_markers[votu]:
                        row += [1]
                        i += 1
                    else:
                        row += [0]
                else:
                    row += [0]
            row += [i]
            master_table[votu]['VOG77'] = i
            out.write( '\t'.join([str(i) for i in row])+'\n')

def write_master_table(master_table):
    '''Write a Master tabel with MGV clustering information '''
    votus_subset = {k:v for k,v in master_table.items() if master_table[k]['MGVref'] != 'MGV' and master_table[k]['decision'] == 'virus'}
    random_votu = list(votus_subset.keys())[0]
    header = list( votus_subset[random_votu].keys() ) 
    with open('master_table_MGV.txt', 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=header, delimiter='\t')
        writer.writeheader()
        for vOTU in votus_subset: writer.writerow(votus_subset[vOTU])


def parse_VOG77(hmmfile):
    '''Parse VOG77 hmmsearch'''
    hits = {}
    VOGs = set()
    with open(hmmfile,'r') as infile:
        for line in infile:
            if line[0]=='#': continue
            r = line.split()
            votu = r[0].rsplit('_',1)[0]
            VOGs.add(r[2])
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
        marker_to_hits[marker_id].append([r[0], r[0].rsplit('_',1)[0] ])
    votu_markers = {}
    for vog,markers in marker_to_hits.items():
        for protein, votu in markers:
            if votu not in votu_markers: votu_markers[votu] = set()
            votu_markers[votu].add(vog)
    return VOGs, votu_markers

def viral_genus_taxonomy(votus_in_genus, master_table):
    '''Pile up host-annotations at different taxa steps , calculate consensus annotation'''
    lineage = ['superkingdom','phylum','class','order','family','genus','species']
    genus_annotation = {'MGV':0,'vOTU':0, 'novel_vOTU':0, 'total_vOTU': 0}
    dominant_host_lineage =  {k:'NA' for k in lineage}
    host_lineage = {k:[] for k in lineage}
    for vOTU in votus_in_genus:
        if vOTU in master_table:   
            if 'MGV' in vOTU: genus_annotation['MGV'] += 1 
            else: 
                if master_table[vOTU]['MGVref'] == 'MGVref': genus_annotation['vOTU'] += 1 
                else: genus_annotation['novel_vOTU'] += 1
                
            if 'hosttax' in master_table[vOTU]:
                vOTU_host_lineage =  master_table[vOTU]['hosttax']
                for i,k in enumerate(lineage):
                    lin = vOTU_host_lineage[i]
                    if lin != 'NA': host_lineage[k].append(lin)
            elif 'host' in master_table[vOTU]: # MGV host information 
                vOTU_host_lineage = master_table[vOTU]['host'].split(';')
                if len(vOTU_host_lineage) == 1: continue
                for i,k in enumerate(lineage):
                    lin = vOTU_host_lineage[i]
                    lin = lin.split('_')[0]
                    if lin != 'NA' and lin != 'NULL': host_lineage[k].append(lin)
    tax_count_string = []
    for k in lineage:
        c = collections.Counter(host_lineage[k])
        cmc = c.most_common()
        if len(cmc) != 0:
            if cmc[0][1]/len(host_lineage[k]) >= 0.25:
                dominant_host_lineage[k] = cmc[0][0]
                taxcount = host_lineage[k].count(cmc[0][0])
                totalannotations = len( host_lineage[k])
                tax_count_string.append( cmc[0][0] + '(' + str(taxcount) + '/' + str(totalannotations) + ')' )
            else:
                tax_count_string.append( cmc[0][0] + '(' + str(taxcount) + '/' + str(totalannotations) + ')' )
    consensus_host_counts = ';'.join(tax_count_string)
    consensus_host_lineage = dominant_host_lineage
    vOTU_crispr_annotated =len(host_lineage['superkingdom'])
    genus_annotation['vOTU_crispr_annotated'] = vOTU_crispr_annotated
    genus_annotation['consensus_host_counts'] = consensus_host_counts 
    genus_annotation['consensus_host_lineage'] = consensus_host_lineage
    genus_annotation['total_vOTU'] = genus_annotation['novel_vOTU'] + genus_annotation['vOTU'] + genus_annotation['MGV'] 
    return genus_annotation


def write_genus_summary_table(viral_genus_summary, fileout='viral_genus_summary.tsv'):
    ''''''
    lineage = ['superkingdom','phylum','class','order','family','genus','species']
    header = ['vOTU','total_vOTU','novel_vOTU','Sardinian_vOTUs','MGV','vOTU_crispr_annotated','consensus_host_counts']
    with open(fileout,'w') as out:
        out.write( '\t'.join( ['viral_genus'] + header + lineage)+'\n')
        for G in viral_genus_summary.keys():
            lineout = [G]
            for i in header:
                lineout += [str(viral_genus_summary[G][i])]
            for lin in lineage:
                lineout += [viral_genus_summary[G]['consensus_host_lineage'][lin]]
            out.write('\t'.join(lineout)+'\n')


### Load Master-tabels for MGV-otus and Dataset specific Viral OTUs.
master_table = tree_functions.load_MGV_master_table(master_table_file= 'master_tabel_final.txt', MGV_file='mgv_contig_info.tsv')
with open('mgv_host_assignments.expanded.tsv','r') as csvfile:
    reader = csv.DictReader(csvfile,delimiter='\t')
    for row in reader:
        MGVID = row['contig_id']
        host = []
        for k in ['superkingdom','phylum','class','order','family','genus','species']:
            host += [row[k]]
        master_table[row['contig_id']]['host'] = ';'.join(host)
master_table = tree_functions.load_vOTU_hosttax(master_table,hosttaxfile='virus.hosttax.txt')


vOTU_enrichment = tree_functions.read_vOTU_enrichment(directory='../../Tables')
VC, virus_VC_lookup = load_species_clusters(master_table,species_file='virus.clusters.95')
viral_genus_clusters = load_genus_clusters(master_table,genus_file='genus_clusters.tsv')
VOGs, votu_markers = parse_VOG77(hmmfile='hmmsearch.txt')

votu_type_counts = {'MGVref':0,'Novel':0}
### 
for vOTU in master_table:
    if master_table[vOTU]['decision'] == 'virus' and master_table[vOTU]['MGVref'] == 'Novel':
        votu_type_counts['Novel'] += 1
    elif master_table[vOTU]['decision'] == 'virus' and master_table[vOTU]['MGVref'] == 'MGVref':
        votu_type_counts['MGVref'] += 1

### Calculte host-consensus across MGV and assembled vOTUs
viral_genus_summary = dict()
for genusid in viral_genus_clusters.keys():
    votus_in_genus = [ v for v in viral_genus_clusters[genusid] if master_table[v]['decision'] == 'virus' ]
    genus_annotation =  viral_genus_taxonomy( viral_genus_clusters[genusid], master_table= master_table)
    viral_genus_summary[genusid] = genus_annotation
    viral_genus_summary[genusid]['Sardinian_vOTUs'] = 0

### Count vOTUs also found in Sardinian Centenarians
for vOTU in master_table.keys():
    if vOTU in vOTU_enrichment:
        S_CE_prevalence = float(vOTU_enrichment[vOTU]['S_CE_prevalence'])
        if S_CE_prevalence > 0:
            genusid = master_table[vOTU]['genus']
            viral_genus_summary[genusid]['Sardinian_vOTUs'] += 1


### Write Viral genus cluster annotation
write_genus_summary_table(viral_genus_summary,fileout='viral_genus_summary.tsv')
 
### Write VOG counts of vOTUs (and MGV representatives)
write_VOG_table(master_table,VOGs,votu_markers)

###
write_master_table(master_table)


### ### ### ###  Tree - Analysis ### ### ### ###  



### Parse tree first time only - attach node names
from Bio import SeqIO
from Bio import Phylo
# tree = Phylo.read("concat.faa.treefile", "newick")
### Name nodes
# n = 0 
# for node in tree.get_nonterminals():
#    n+=1
#    node.name = 'n'+str(n)
# Phylo.write(tree,"new.tree.nwk","newick")

tree = Phylo.read("new.tree.nwk", "newick")
leafs = tree.get_terminals()

master_table = tree_functions.load_MGV_master_table(master_table_file= 'master_table_MGV.txt', MGV_file='mgv_contig_info.tsv')
master_table = tree_functions.load_vOTU_hosttax(master_table,hosttaxfile='virus.hosttax.txt')

black_listed = [l.name.replace('__','::') for l in leafs if not l.name.replace('__','::') in master_table]
leafs = [c for c in leafs if not c.name.replace('__','::') in black_listed ]


### 
viral_genus_dict = {'MGV':set(),'Novel':set()}
tree_counts = {'MGVref':0,'Novel':0}
for l in leafs:
    if not 'MGV' in l.name:
        vOTU = l.name.replace('__','::')
        if master_table[vOTU]['decision'] == 'virus' and master_table[vOTU]['MGVref'] == 'Novel':
            tree_counts['Novel'] += 1
            viral_genus_dict['Novel'].add(master_table[vOTU]['genus'])
        elif master_table[vOTU]['decision'] == 'virus' and master_table[vOTU]['MGVref'] == 'MGVref':
            tree_counts['MGVref'] += 1
            viral_genus_dict['MGV'].add(master_table[vOTU]['genus'])

actual_novel_genera = viral_genus_dict['Novel'].difference(viral_genus_dict['MGV'])
shared_MGV_genera = viral_genus_dict['Novel'].intersection(viral_genus_dict['MGV'])

print(len(actual_novel_genera),' new genus of virus in Tree')
print(len(shared_MGV_genera),' genus shared with MGV of virus in Tree')

### Relevant research questions 
# Can we place any Viruses that actually form new Genus? 
# In what Genera (with MGV references) are the rest placed - colour branches by MCRA

clades = tree_functions.lookup_by_names(tree)
leafs_by_genus = {x:[] for x in list(viral_genus_dict['MGV'])+list(viral_genus_dict['Novel']) }
for leaf in leafs:
    if not 'MGV' in leaf.name: 
        vOTU = leaf.name.replace('__','::')
        if not vOTU in black_listed:
            genus = master_table[vOTU]['genus']
            leafs_by_genus[genus]+=[leaf]
    else:
        genus = master_table[vOTU]['genus']
        if genus in leafs_by_genus: leafs_by_genus[genus]+=[leaf]

genus_to_mcra_nodes, genus_prevalence = tree_functions.genus_annotation(tree,leafs_by_genus,directory='../../Tables')
genus_abundance = tree_functions.read_genus_abundance(directory='../../Tables')
vOTU_prevalence = tree_functions.read_vOTU_prevalence(directory='../../Tables')
vOTU_enrichment = tree_functions.read_vOTU_enrichment(directory='../../Tables')



### Summarise by vOTUs by genera  
for vOTU in master_table.keys():
    if not vOTU in vOTU_enrichment: continue
    ### Enrichment 
    CE_young, CE_elderly, elderly_young = vOTU_enrichment[vOTU]['CE_young_adjusted'], vOTU_enrichment[vOTU]['CE_elderly_adjusted'], vOTU_enrichment[vOTU]['elderly_young_adjusted']

    ### Sardinian Prevalence 
    S_CE_prevalence = float(vOTU_enrichment[vOTU]['S_CE_prevalence'])
    master_table[vOTU]['S_CE_prevalence'] = S_CE_prevalence
    if float(CE_young) < 0.05 and float(CE_elderly) < 0.05:
        master_table[vOTU]['CE_enriched'] = True
    else:
        master_table[vOTU]['CE_enriched'] = False



### ### ### ### Tree annotation ### ### ### ### 
os.makedirs('tree_annotation',exist_ok = True)
tree_functions.prune(leafs,directory='tree_annotation')

tree_functions.write_hostphylum_color(leafs,master_table,viral_genus_summary,directory='tree_annotation')
tree_functions.write_novely_indicator(leafs,master_table, directory='tree_annotation')
tree_functions.write_vOTU_prevalence_bar(leafs,vOTU_prevalence,master_table,directory='tree_annotation')
tree_functions.write_vOTU_enrichment_tiles(leafs,vOTU_enrichment,directory='tree_annotation')
tree_functions.write_CEgenus_color(leafs,master_table,directory='tree_annotation')








