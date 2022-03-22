#!/bin/python
import os
import csv
def load_vOTU_hosttax(master_table,hosttaxfile):
    tax = ['superkingdom','phylum','class','order','family','genus','species']
    with open(hosttaxfile) as csvfile:
        reader = csv.DictReader(csvfile,delimiter='\t')
        for row in reader:
            if not row['virus'] in master_table: continue 
            taxonomy = []
            for taxlevel in tax:
                taxonomy += [row[taxlevel]]
            master_table[row['virus']]['hosttax'] = taxonomy
    return master_table

def load_MGV_master_table(master_table_file, MGV_file):
    master_table = {}
    with open(master_table_file) as csvfile:
        reader = csv.DictReader(csvfile,delimiter='\t')
        for row in reader:
            master_table[row['contig_id']] = row
    mgv_annotation = {}
    with open(MGV_file) as csvfile:
        reader = csv.DictReader(csvfile,delimiter='\t')
        for row in reader:
            row['MGVref'] = 'MGV'
            row['genus'] = 'NA'
            row['decision'] = 'virus'
            row['host'] = 'NA'
            mgv_annotation[row['contig_id']] = row
    master_table.update(mgv_annotation)
    return master_table


def read_genus_abundance(directory):
    file = os.path.join(directory,'Virus_mean_TPM.tsv')
    genus_abundance = dict()
    with open(file) as csvfile:
        reader = csv.DictReader(csvfile,delimiter='\t')
        for row in reader:
            row['Centenarian'] = round(float(row['Centenarian']),2)
            row['Elderly_control'] = round(float(row['Elderly_control']),2)
            row['Young_control'] = round(float(row['Young_control']),2)
            genus_abundance[row['genus']] = row
    return genus_abundance

def read_vOTU_prevalence(directory):
    file = os.path.join(directory,'Virus_prevalence.tsv')
    vOTU_prevalence = dict()
    with open(file) as csvfile:
        reader = csv.DictReader(csvfile,delimiter='\t')
        for row in reader:
            row['Centenarian'] = round(float(row['Centenarian']),2)
            row['Elderly_control'] = round(float(row['Elderly_control']),2)
            row['Young_control'] = round(float(row['Young_control']),2)
            vOTU_prevalence[row['virus']] = row
    return vOTU_prevalence


def read_vOTU_enrichment(directory):
    file = os.path.join(directory,'Virus_enrichment.tsv')
    vOTU_enrichment = dict()
    with open(file) as csvfile:
        reader = csv.DictReader(csvfile,delimiter='\t')
        for row in reader:
            vOTU_enrichment[row['virus']] = row
    return vOTU_enrichment



def genus_annotation(tree,leafs_by_genus,directory):
    file = os.path.join(directory,'Virus_genus_prevalence.tsv')
    genus_prevalence = dict()
    with open(file) as csvfile:
        reader = csv.DictReader(csvfile,delimiter='\t')
        for row in reader:
            row['Centenarian'] = round(float(row['Centenarian']),2)
            row['Elderly_control'] = round(float(row['Elderly_control']),2)
            row['Young_control'] = round(float(row['Young_control']),2)
            genus_prevalence[row['genus']] = row
    
    ### Get MCRA clade for each genus
    if not os.path.exists('genus_mcra_node.txt'):
        print('Calculating genus-wise MCRA nodes')
        genus_to_mcra = dict()
        for genus,leafs in leafs_by_genus.items():
            if len(leafs) == 1:
                genus_to_mcra[genus] = get_parent(tree,leafs[0])
            else:
                genus_to_mcra[genus] = mcra(tree,leafs)
        genus_to_mcra_nodes = {genus:clade.name for genus,clade in genus_to_mcra.items()}
        with open('genus_mcra_node.txt','w') as outfile:
            for genus, node in genus_to_mcra_nodes.items():
                outfile.write('{}\t{}\n'.format(genus,node))
    
    genus_to_mcra_nodes = {}
    with open('genus_mcra_node.txt','r') as filein:
        for line in filein:
            genus,node = line.strip().split()
            genus_to_mcra_nodes[genus] = node
    return genus_to_mcra_nodes,genus_prevalence

def lookup_by_names(tree):
    names = {}
    for clade in tree.find_clades():
        if clade.name:
            if clade.name in names:
                raise ValueError("Duplicate key: %s" % clade.name)
            names[clade.name] = clade
    return names

def get_parent(tree, child_clade):
    '''Get parent of child clade '''
    node_path = tree.get_path(child_clade)
    return node_path[-2]

def mcra(tree,leafs):
    mcra = tree.common_ancestor(leafs)
    return mcra






############  Tree annotation functions ############
def prune(leafs,directory):

    fileout= os.path.join(directory,'tree_annotation.prune.txt')
    out = open(fileout,'w')
    out.write('PRUNE\nDATA\n')
    non_ref= []
    for leaf in leafs:
        if not 'MGV' in leaf.name: 
            out.write('{}\n'.format(leaf.name))
    out.close()


def write_hostphylum_color(leafs,master_table, viral_genus_summary, directory):

    fileout= os.path.join(directory,'tree_annotation.hostphylum.txt')
    out = open(fileout,'w')
    out.write('DATASET_COLORSTRIP\nSEPARATOR TAB\n')
    out.write('DATASET_LABEL\tHostphylum\n')
    out.write('STRIP_WIDTH\t25\n')
    out.write('DATA\n')

    non_ref= []
    for leaf in leafs:
        if not 'MGV' in leaf.name: non_ref.append(leaf)
    
    phylum_colors = {'Firmicutes':'#009E73','Actinobacteria':'#D55E00','Actinobacteriota':'#D55E00',
    'Bacteroidetes':'#F0E442','Bacteroidota':'#F0E442',
    'Proteobacteria':'#0072B2',
    'Verrucomicrobia':'#CC79A7',
    'NA':'#f7f7f5',
    'Desulfobacterota':'#92175B'}
    for leaf in non_ref:
        actual_name = leaf.name.replace('__','::')
        virus_genus = master_table[actual_name]['genus']
        phylum = viral_genus_summary[virus_genus]['consensus_host_lineage']['phylum']
        #phylum = master_table[actual_name]['hosttax'][1]
        if not phylum in phylum_colors:
            print(phylum)
            phylum = 'NA'
        out.write('{}\t{}\n'.format(leaf.name,phylum_colors[phylum]))
    out.close()



def write_CEgenus_color(leafs,master_table, directory):

    fileout= os.path.join(directory,'tree_annotation.viralgenus.txt')
    out = open(fileout,'w')
    out.write('DATASET_COLORSTRIP\nSEPARATOR TAB\n')
    out.write('DATASET_LABEL\tviralgenus\n')
    out.write('STRIP_WIDTH\t25\n')
    out.write('DATA\n')

    non_ref= []
    for leaf in leafs:
        if not 'MGV' in leaf.name: non_ref.append(leaf)
    
    genus_colors = {'G7':'#E69F00','G78':'#AD8E49','G2':'#FF5852','NA':'#f7f7f5'}
    for leaf in non_ref:
        actual_name = leaf.name.replace('__','::')
        genus = master_table[actual_name]['genus']
        if not genus in genus_colors:
            genus = 'NA'
        out.write('{}\t{}\n'.format(leaf.name,genus_colors[genus]))
    out.close()



def write_dataset_indicator(leafs,master_table,directory):

    fileout= os.path.join(directory,'tree_annotation.dataset.txt')
    out = open(fileout,'w')
    out.write('DATASET_COLORSTRIP\nSEPARATOR TAB\n')
    out.write('DATASET_LABEL\tDataset\n')
    out.write('DATA\n')
    out.write('CONTAINS==MGV\t#080707\n')
    non_ref= []
    for leaf in leafs:
        if not 'MGV' in leaf.name: non_ref.append(leaf)
    
    for leaf in non_ref:
        actual_name = leaf.name.replace('__','::')
        out.write('{}\t#2c911a\n'.format(leaf.name))
    out.close()


def write_novely_indicator(leafs,master_table, directory):
    fileout= os.path.join(directory,'tree_annotation.novel.txt')
    out = open(fileout,'w')
    out.write('DATASET_COLORSTRIP\nSEPARATOR TAB\n')
    out.write('DATASET_LABEL\tNovel\n')
    out.write('COLOR_BRANCHES\t1\n')
    out.write('DATA\n')

    non_ref= []
    for leaf in leafs:
        if not 'MGV' in leaf.name: non_ref.append(leaf)
    
    for leaf in non_ref:
        actual_name = leaf.name.replace('__','::')
        if master_table[actual_name]['MGVref'] == 'Novel': out.write('{}\t#c40210\n'.format(leaf.name))
        else: out.write('{}\t#020d01\n'.format(leaf.name))
    out.close()



def write_vOTU_enrichment_tiles(leafs,vOTU_enrichment,directory):
    fileout= os.path.join(directory,'tree_annotation.vOTUenrichment.txt')
    out = open(fileout,'w')
    out.write('DATASET_HEATMAP\nSEPARATOR TAB\n')
    out.write('DATASET_LABEL\tWilcox_Enrichment\n')
    out.write('COLOR\t#ff0000\n')
    out.write('FIELD_SHAPES\t2\t2\t2\n')
    out.write('FIELD_LABELS\tCE_young\tCE_elderly\telderly_young\n')
    out.write('STRIP_WIDTH\t25\n')
    out.write('AUTO_LEGEND\t0\n')
    out.write('SHOW_TREE\t0\n')
    out.write('BORDER_WIDTH\t1\n')
    out.write('BORDER_COLOR\t#eeeeee\n')
    out.write('COLOR_MIN\t#ffffff\n')
    out.write('COLOR_MAX\t#0a0a0a\n')
    out.write('DATA\n')
    non_ref= []
    for leaf in leafs:
        if not 'MGV' in leaf.name: non_ref.append(leaf)
    for leaf in non_ref:
        actual_name = leaf.name.replace('__','::')
        CE_young, CE_elderly, elderly_young = vOTU_enrichment[actual_name]['CE_young_adjusted'], vOTU_enrichment[actual_name]['CE_elderly_adjusted'], vOTU_enrichment[actual_name]['elderly_young_adjusted']
        if float(CE_young) <0.05: CE_young = 1
        else: CE_young =0
        if float(CE_elderly) <0.05: CE_elderly = 1 
        else: CE_elderly =0
        if float(elderly_young) <0.05: elderly_young = 1 
        else: elderly_young =0
        out.write('{}\t{}\t{}\t{}\n'.format(leaf.name,CE_young,CE_elderly,elderly_young))
    out.close()

def write_vOTU_prevalence_bar(leafs,vOTU_prevalence,master_table, directory):

    fileout= os.path.join(directory,'tree_annotation.vOTUprevalence.txt')
    out = open(fileout,'w')
    out.write('DATASET_MULTIBAR\nSEPARATOR TAB\n')
    out.write('DATASET_LABEL\tAge_prevalence\n')
    out.write('FIELD_COLORS\t#f0a207\t#1d77ad\t#bfc4c7\n') #Centenarian, Elderly, Young that's the order
    out.write('FIELD_LABELS\tCE\tElderly\tYoung\n')
    #out.write('DATASET_SCALE\t1\t1\t1\n')
    out.write('ALIGN_FIELDS\t1\n')
    out.write('WIDTH\t200\n')
    out.write('DATA\n')
    # ID1,value1,value2,value3
    non_ref= []
    for leaf in leafs:
        if not 'MGV' in leaf.name: non_ref.append(leaf)
    
    for leaf in non_ref:
        actual_name = leaf.name.replace('__','::')
        CE, elderly, young = vOTU_prevalence[actual_name]['Centenarian'], vOTU_prevalence[actual_name]['Elderly_control'], vOTU_prevalence[actual_name]['Young_control']
        out.write('{}\t{}\t{}\t{}\n'.format(leaf.name,CE,elderly,young))
    out.close()






def write_genus_prevalence_bar(tree,leafs,genus_to_mcra_nodes,genus_prevalence,master_table, directory):

    fileout= os.path.join(directory,'tree_annotation.genusprevalence.txt')
    out = open(fileout,'w')
    out.write('DATASET_MULTIBAR\nSEPARATOR TAB\n')
    out.write('DATASET_LABEL\tAge_prevalence\n')
    out.write('FIELD_COLORS\t#f0a207\t#1d77ad\t#bfc4c7\n') #Centenarian, Elderly, Young that's the order
    out.write('FIELD_LABELS\tCE\tElderly\tYoung\n')
    out.write('DATASET_SCALE\t1\t1\t1\n')
    out.write('ALIGN_FIELDS\t1\n')
    out.write('WIDTH\t100\n')
    out.write('DATA\n')
    # ID1,value1,value2,value3

    non_ref= []
    for leaf in leafs:
        if not 'MGV' in leaf.name: non_ref.append(leaf)
        for leaf in non_ref:
            actual_name = leaf.name.replace('__','::')
            genus = master_table[actual_name]['genus']
            genus_node = genus_to_mcra_nodes[genus]
            novelty = genus_prevalence[genus]['MGVref']
    
    for leaf in non_ref:
        actual_name = leaf.name.replace('__','::')
        genus = master_table[actual_name]['genus']
        CE, elderly, young = genus_prevalence[genus]['Centenarian'], genus_prevalence[genus]['Elderly_control'], genus_prevalence[genus]['Young_control']
        out.write('{}\t{}\t{}\t{}\n'.format(leaf.name,CE,elderly,young))
    out.close()


def write_viralfamily_color(leafs,master_table, directory):

    fileout= os.path.join(directory,'tree_annotation.viralfamily.txt')
    out = open(fileout,'w')
    out.write('DATASET_COLORSTRIP\nSEPARATOR TAB\n')
    out.write('DATASET_LABEL\tViralfamily\n')
    out.write('STRIP_WIDTH\t25\n')
    out.write('DATA\n')
    non_ref= []
    for leaf in leafs:
        if not 'MGV' in leaf.name: non_ref.append(leaf)
    
    phylum_colors = {'crAss-phage':'#a80398','Siphoviridae':'#f5bd3b','Myoviridae':'#3be9f5','Microviridae':'#c40613','Podoviridae':'#420307','NA':'#f7f7f5'}
    for leaf in non_ref:
        actual_name = leaf.name.replace('__','::')
        phylum = master_table[actual_name]['ictv_family']
        if not phylum in phylum_colors:
            phylum = 'NA'
        out.write('{}\t{}\n'.format(leaf.name,phylum_colors[phylum]))
    out.close()
