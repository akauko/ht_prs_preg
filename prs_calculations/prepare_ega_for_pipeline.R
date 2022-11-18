#!/usr/bin/env Rscript

#args = commandArgs(trailingOnly=TRUE)

## test if there is two argument: if not, return an error
#if (length(args)!=3) {
#  stop("Two arguments must be supplied (input file).n", call.=FALSE)
#} 

#file_in <- args[1]
#file_rsin <- args[2]
#file_out <- args[3]

data_path <- "~/ht_prs_preg/data/ukb_ega_gwas/"

#file_in <- paste0(data_path, "test.tsv")
file_in <- paste0(data_path, "mat_all_chrALL_STERR_EU_exclFIN_VS_ASIA.1tbl.gz")
file_posin <- paste0(data_path, "meta_variants_chr_pos")
file_out <- paste0(data_path, "ega2.pe.gwas.tsv.gz")


library(data.table, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
library(magrittr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(stringr, warn.conflicts = FALSE)

#Changes chromosome positions to rsid's, extracts columns for alleles and removes unnecessary columns
#This processes only one file. Looping at bash. These are large files and processing even one takes time.

### position mapping is read in
posin <- fread(file_posin, header = F) %>%
  rename(MarkerName = V1, Chromosome = V2, Position = V3)


### gwas summary
df <- fread(file_in, header = T) %>% 
  left_join(posin, by="MarkerName") %>%
  mutate(rsid = str_extract(MarkerName, "rs\\d+")) %>%
  mutate(MarkerName = if_else(is.na(rsid), MarkerName, rsid)) %>%
  mutate_at(c("Allele1","Allele2"), str_to_upper)
  #The file did not have X-chromosome, so no need to remove it.
  #NA's present only in rsid mapping, not in other variables, so no NA removal step

#df %>% summarise_all(~sum(is.na(.)))  

## New file is written out
fwrite(df, file_out, sep="\t")

