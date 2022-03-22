#!/usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE)
print(args)
library(tidyverse)
library(dplyr)
library(readr)
library(data.table)
library(fastmatch)
### We need Fastmatch 
#install.packages('fastmatch')
#require('fastmatch')

### Arguments 
# The programme expects files in the following order like
# genecat.virus.blast.m6 genecat.TPM.matrix.tsv directoryout 
### 
outdir = args[3]
TPMfile = args[2]

### Filter for genes with at least 95% sequence identity
blast_table <- fread(args[1]) %>% as_tibble() 
blast_table <- blast_table[blast_table$V3>=95,]

### Listify genes belonging to each Viralbin 
virgenes <- lapply(unique(blast_table$V2),function(vir){
  tmp <- blast_table[blast_table$V2==vir,]
  tmp <- tmp[,c('V1','V2')]
  return(tmp)
})
virtable <- do.call(rbind,virgenes) %>% select(V2,V1) 


vir_gene_list <- tapply(1:nrow(virtable), virtable$V2, function(x){
  return(as.vector(virtable$V1[x]))
})
vir_gene_list <- vir_gene_list[sapply(vir_gene_list, length) > 3]
viruses <- names(vir_gene_list)
# Prune MSP table to only include bins with > 3 genes
virtable <- virtable %>% filter(V2 %in% names(vir_gene_list))


### Load TPM matrix
sample_tpm <- read_delim(TPMfile, delim ='\t', col_names=TRUE, progress  = TRUE)
samples  <- colnames(sample_tpm)
rownames <- unlist(sample_tpm[,1])
sample_tpm_subset <- sample_tpm[rownames %in% blast_table$V1,]
write_tsv(file=str_c(outdir,'/',"TPM_virus_subset.tsv"),sample_tpm_subset, col_names=TRUE)
rm(sample_tpm_subset)

cat("Calculating medians for each Virus \n")
TPM_sample_medians <- lapply(samples[-1], function(i){
  cat(i,'\n')
  sample_index <- fmatch(i, names(sample_tpm))
  sample_tpm_column <- sample_tpm[,c(1,sample_index)]

  TPM_median <- lapply(vir_gene_list, function(x){
    indices <- fmatch(x,sample_tpm_column[[1]])
    vals <- sample_tpm_column[[2]][indices]
    medianAbundance <- median(vals)
    return(medianAbundance)
  })
  TPM_median_df <- plyr::ldply(TPM_median)
  return(TPM_median_df[[2]])
})
TPM_sample_matrix <- do.call(cbind,TPM_sample_medians)
rownames(TPM_sample_matrix) <- viruses
colnames(TPM_sample_matrix) <- samples[-1]
write_tsv(as.data.frame(cbind(viruses,TPM_sample_matrix)), file=str_c(outdir,'/',"combined_viralTPM.txt"))

rm(TPM_sample_matrix,TPM_sample_medians)


cat("Calculating gene detection frequencies for each Virus \n")
TPM_sample_frequencies <- lapply(samples[-1], function(i){
  cat(i,'\n')
  sample_index <- fmatch(i, names(sample_tpm))
  sample_tpm_column <- sample_tpm[,c(1,sample_index)]

  TPM_freq <- lapply(vir_gene_list, function(x){
  indices <- fmatch(x,sample_tpm_column[[1]])
  vals <- sample_tpm_column[[2]][indices]
  freq <- sum(vals>0)/length(vals)
  })
  TPM_freq_df <- plyr::ldply(TPM_freq)
  return(TPM_freq_df[[2]])
})
Freq_sample_matrix <- do.call(cbind,TPM_sample_frequencies)
rownames(Freq_sample_matrix) <- viruses
colnames(Freq_sample_matrix) <- samples[-1]

write_tsv(as.data.frame(cbind(viruses,Freq_sample_matrix)), file=str_c(outdir,'/',"combined_viralFrequency.txt"))
write_delim(virtable , file= str_c(outdir,'/','virus_gene_table.tsv'),delim='\t')
cat("Done writing results. Exiting R\n")

