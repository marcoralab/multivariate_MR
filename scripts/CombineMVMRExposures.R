library(tidyr)
library(readr)
library(dplyr)

exp1_mrdat.path = snakemake@input[["exp1_mrdat"]] # Harmonized MR data
exp2_mrdat.path = snakemake@input[["exp2_mrdat"]] # Output
out.path = snakemake@output[["out"]] # Output

exp1_mrdat <- read_csv(exp1_mrdat.path)
exp2_mrdat <- read_csv(exp2_mrdat.path)

bind_rows(exp1_mrdat, exp2_mrdat) %>%
  filter(mr_keep == TRUE, Outliers == FALSE) %>%
  distinct(SNP) %>%
  select(SNP) %>%
  write_tsv(out.path, col_names = T)
