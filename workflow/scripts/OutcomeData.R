#!/usr/bin/Rscript

### ===== Command Line Arguments ===== ##
exposure.summary = snakemake@input[["ExposureSummary"]] # Exposure summary statistics
outcome.summary = snakemake@input[["OutcomeSummary"]] # Outcome Summary statistics


### ===== Load packages ===== ###
library(tidyverse)   ## For data wrangling
#suppressMessages(library(Hmisc))       ## Contains miscillaneous funtions

### ===== READ IN SNPs ===== ###
message("READING IN EXPOSURE \n")
exposure.dat <- read_tsv(exposure.summary)

message("\n READING IN OUTCOME \n")
outcome.dat.raw <- read_tsv(outcome.summary, guess_max = 15000000)

### ===== EXTACT SNPS ===== ###
message("\n EXTRACTING SNP EFFECTS FROM OUTCOME GWAS  \n")
outcome.dat <- outcome.dat.raw %>%
  right_join(select(exposure.dat, SNP), by = 'SNP')


### ===== MISSING SNPS SNPS ===== ###

outcome.dat %>%
  filter(is.na(CHROM)) %>%
  select(SNP) %>%
  write_tsv(snakemake@output[["missing"]], col_names = F)

message("\n EXPORTING \n")
## Write out outcomes SNPs
write_tsv(outcome.dat, snakemake@output[["snps"]])
