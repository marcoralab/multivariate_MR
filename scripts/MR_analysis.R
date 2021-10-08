#!/usr/bin/Rscript
## Run MR analysis and Senstivity analysis for w/ outliers retained and removed
if(any(grepl("conda", .libPaths(), fixed = TRUE))){
  message("Setting libPaths")
  df = .libPaths()
  conda_i = which(grepl("conda", df, fixed = TRUE))
  .libPaths(c(df[conda_i], df[-conda_i]))
}
.libPaths(c(.libPaths(), "/hpc/users/harern01/.Rlib"))

library(tidyr)
library(readr)
library(dplyr)
library(forcats)
library(TwoSampleMR) ## For conducting MR https://mrcieu.github.io/TwoSampleMR/

input = snakemake@input[["mrdat"]] # Harmonized MR data
output = snakemake@params[["out"]] # Output

mrdat.raw <- read_csv(input)
pt = mrdat.raw %>% slice(1) %>% pull(pt)
n_ambiguous = nrow(mrdat.raw) - nrow(mrdat.raw %>% filter(mr_keep == TRUE))
n_gws <- nrow(mrdat.raw %>% filter(mr_keep == TRUE)) - nrow(mrdat.raw %>% filter(mr_keep == TRUE, pval.outcome > 5e-8))
n_outliers <- nrow(mrdat.raw %>% filter(mr_keep == TRUE, pval.outcome > 5e-8)) - nrow(mrdat.raw %>% filter(mr_keep == TRUE, pval.outcome > 5e-8, Outliers == FALSE))

## Remove SNPs that are:
#   palindromic
#   ambigous
#   pval.outcome < 5e-8
#   outlier
mrdat <- mrdat.raw %>%
  filter(mr_keep == TRUE, pval.outcome > 5e-8, Outliers == FALSE)

nsnps <- tibble(
  id.exposure = as.character(mrdat.raw[1,'id.exposure']),
  id.outcome = as.character(mrdat.raw[1,'id.outcome']),
  outcome = as.character(mrdat.raw[1,'outcome']),
  exposure = as.character(mrdat.raw[1,'exposure']),
  pt = pt,
  nsnps = nrow(mrdat.raw),
  n_ambiguous = n_ambiguous,
  n_outcome_gws = n_gws,
  n_outliers = n_outliers,
  n_mrsnps = nrow(mrdat)
)

## ================= w/o outliers ================= ##
## MR analysis
mr_res <- mr(mrdat, method_list = c("mr_ivw_fe", "mr_weighted_median", "mr_weighted_mode", "mr_egger_regression")) %>%
  as_tibble() %>%
  mutate(pt = pt) %>%
  select(id.exposure, id.outcome, outcome, exposure, pt, method, nsnp, b, se, pval)

## Cochrans Q heterogeneity test
mr_heterogenity <- mr_heterogeneity(mrdat, method_list=c("mr_egger_regression", "mr_ivw")) %>%
  as_tibble() %>%
  mutate(pt = pt) %>%
  select(id.exposure, id.outcome, outcome, exposure, pt, method, Q, Q_df, Q_pval)

## MR Egger Test of Pliotropy
mr_plei <- mr_pleiotropy_test(mrdat) %>%
  as_tibble() %>%
  mutate(pt = pt) %>%
  select(id.exposure, id.outcome, outcome, exposure, pt, egger_intercept, se, pval)

## ================= Write Out ================= ##
mr_heterogenity %>% write_tsv(paste0(output, '_MR_heterogenity.txt'))
mr_plei %>% write_tsv(paste0(output, '_MR_egger_plei.txt'))
mr_res %>% write_tsv(paste0(output, '_MR_Results.txt'))
nsnps %>% write_tsv(paste0(output, '_MR_nsnps.txt'))
