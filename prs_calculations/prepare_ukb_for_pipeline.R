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

#file_in <- paste0(data_path, "test2.tsv")
file_in <- paste0(data_path, "4080_irnt.gwas.imputed_v3.female.tsv.gz")
#file_in <- paste0(data_path, "20002_1073.gwas.imputed_v3.female.tsv.gz")
file_rsin <- paste0(data_path, "variants.tsv.gz")
file_out <- paste0(data_path, "ukb.sbp.female.gwas.tsv.gz")
#file_out <- paste0(data_path, "ukb.20002_1073.female.gwas.tsv.gz")

library(data.table, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
library(magrittr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(stringr, warn.conflicts = FALSE)

#Changes chromosome positions to rsid's, extracts columns for alleles and removes unnecessary columns
#This processes only one file. Looping at bash. These are large files and processing even one takes time.

### variant - rs mapping is read in. Variants without rsid are dropped
rsin <- fread(file_rsin) %>% 
  select(variant, rsid) %>%
  filter(str_detect(rsid, "rs"))

### gwas summary
df <- fread(file_in, header = T, sep = "\t") %>% 
  separate(variant, c("chr", "pos", "ref", "alt"), sep = ":", remove="F") %>%
  left_join(rsin, by="variant") %>%
  mutate(variant = if_else(is.na(rsid), variant, rsid)) %>%
  drop_na(beta, se, pval) %>%
  filter(!str_detect(variant, "^X"))

#df %>% summarise_all(~sum(is.na(.)))  
  
## New file is written out
fwrite(df, file_out, sep="\t")

