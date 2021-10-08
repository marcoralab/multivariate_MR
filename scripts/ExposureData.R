#!/usr/bin/Rscript

### ===== Command Line Arguments ===== ##
exposure.summary = snakemake@input[["summary"]] # Exposure summary statistics
p.threshold = as.numeric(snakemake@params[["Pthreshold"]])
exposure.clump = snakemake@input[["ExposureClump"]]
out.file = snakemake@output[["out"]] # SPECIFY THE OUTPUT FILE

### ===== Load packages ===== ###
library(tidyverse)   ## For data wrangling

### ===== Read in Data ===== ###
message("\n READING IN EXPOSURE \n")
exposure.dat <- read_tsv(exposure.summary, guess_max = 15000000) %>%
  filter(P < p.threshold)

### ===== Clump Exposure ===== ###
message("\n CLUMPING EXPOSURE SNPS \n")

## Plink Pre-clumped
mr_exposure.dat_ld <- read_table2(exposure.clump) %>%
  filter(!is.na(CHR)) %>%
  select(CHR, F, SNP, BP, P,TOTAL,NSIG)

# Filter exposure data for clumped SNPs
exposure.dat <- exposure.dat %>% filter(SNP %in% mr_exposure.dat_ld$SNP)

### ===== Write Out Exposure ===== ###
message("\n Writing Out Exposure \n")

write_tsv(exposure.dat, out.file)
