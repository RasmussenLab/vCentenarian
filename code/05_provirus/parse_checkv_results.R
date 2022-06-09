#!/usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE)

library(tidyverse)
library(data.table)

### Input arguments
quality_file = args[1]
completeness_file = args[2]
outfile_basename = args[3]


checkv_quality <- read_tsv(quality_file,col_names = T)
checkv_completeness <- read_tsv(completeness_file,col_names = T)



### Select High confidence virus bins ###
HQ_bins <- checkv_quality %>% filter(checkv_quality=='High-quality')

HQ_bins_high_confidence <- checkv_completeness %>% filter(contig_id %in% HQ_bins$contig_id) %>% filter(aai_id>75,aai_af>30) %>%
  mutate(length_difference_percentage= (abs(aai_expected_length-contig_length)/aai_expected_length)*100 ) %>%
  filter(length_difference_percentage<15) %>% select(-length_difference_percentage) %>% pull(contig_id) 

### Add Circular bins as well
circular_bins <- checkv_quality %>% filter(checkv_quality=='Complete',completeness_method=='DTR (high-confidence)')  
circular_bins_high_confidence <- checkv_completeness %>% filter(contig_id %in% circular_bins$contig_id)  %>% filter(aai_confidence=='high') %>% pull(contig_id)

HQ_high_confidence_viruses <- rbind(checkv_quality %>% filter(contig_id %in% HQ_bins_high_confidence),  checkv_quality %>% filter(contig_id %in% circular_bins_high_confidence )  )

fileout = paste0(outfile_basename,'.HQ','.tsv')
HQ_high_confidence_viruses  %>% mutate(viraltype = 'HQ') %>% write_tsv(fileout,col_names = T)


### Select also Medium quality high confidence virus bins ###
MQ_bins <- checkv_quality %>% filter(checkv_quality=='Medium-quality')

MQ_bins_high_confidence <- checkv_completeness %>% filter(contig_id %in% MQ_bins$contig_id) %>% filter(aai_id>75,aai_af>30) %>%
  mutate(length_difference_percentage= (abs(aai_expected_length-contig_length)/aai_expected_length)*100 ) %>%
  filter(length_difference_percentage<15) %>% select(-length_difference_percentage) %>% pull(contig_id) 

fileout = paste0(outfile_basename,'.MQ','.tsv')

checkv_quality %>% filter(contig_id %in% MQ_bins_high_confidence) %>% mutate(viraltype = 'MQ') %>% write_tsv(fileout,col_names= T)

### Select Confident Pro-virus predictions

Provirus_bins_high_confidence <- checkv_quality %>% filter(provirus == 'Yes', !grepl('HMM-based',completeness_method), checkv_quality %in% c('Medium-quality','High-quality','Complete')) %>% pull(contig_id)
fileout = paste0(outfile_basename,'.provirus','.tsv')
checkv_quality %>% filter(contig_id %in% Provirus_bins_high_confidence) %>% mutate(viraltype = 'provirus') %>% write_tsv(fileout,col_names= T)

### More dark-matter like bins , little host evidence
LQ_darkmatter <- checkv_quality %>% filter(completeness<50) %>% arrange(-viral_genes) %>% filter(host_genes<3, viral_genes>=5,viral_genes>host_genes,gene_count>=10)

fileout = paste0(outfile_basename,'.LQ','.tsv')

LQ_darkmatter %>%  mutate(viraltype = 'LQ') %>%  write_tsv(fileout,col_names= T)







