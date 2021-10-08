#!/usr/bin/Rscript
if(any(grepl("conda", .libPaths(), fixed = TRUE))){
  message("Setting libPaths")
  df = .libPaths()
  conda_i = which(grepl("conda", df, fixed = TRUE))
  .libPaths(c(df[conda_i], df[-conda_i]))
}
.libPaths(c(.libPaths(), "/hpc/users/harern01/.Rlib"))
### ===== Command Line Arguments ===== ##
infile = snakemake@input[["mrdat"]] # Exposure summary statistics
out = snakemake@params[["out"]]

### ===== Load packages ===== ###
library(tidyr)
library(readr)
library(dplyr)
library(forcats)
library(TwoSampleMR) ## For formating data
library(RadialMR) ## For detecting pleitropy

message('\nREAD IN DATA\n')
mrdat.raw <- read_csv(infile)

## Format data for RadialMR
## Run IVW Radial Regression
message('\nFORMAT DATA FOR RADIAL MR\n')
radial_dat <- (
  mrdat.raw %>%
    filter(pval.outcome > 5e-8) %>%
    filter(mr_keep == TRUE) %>%
    dat_to_RadialMR(.)
)[[1]]

message('\nRUNNING RADIAL MR\n')
sink(paste0(out, '_radialMRoutput.txt'), append=FALSE, split=FALSE)
radial.out <- ivw_radial(radial_dat, weights = 3, alpha = 0.05/nrow(radial_dat))
sink()

## Merge Outliers with MR Data
mrdat.out <- mrdat.raw %>%
  left_join(select(radial.out$data, SNP, Qj, Qj_Chi, Outliers), by = 'SNP') %>%
  rename(Q_statistic = Qj, Q_pval = Qj_Chi) %>%
  mutate(Outliers = fct_recode(Outliers, 'FALSE' = 'Variant', 'TRUE' = 'Outlier'), Outliers = as.logical(Outliers))

## Extracting Radial regression MR estimates
# radial.ivw <- as_tibble(radial.out$coef) %>%
#   mutate(model = 'Effect (Mod.2nd)') %>%
#   rename(b = Estimate, se = 'Std. Error', t.stat = 't value', p = 'Pr(>|t|)')
# radial.ivw.iterative <- as_tibble(radial.out$it.coef) %>%
#   mutate(model = 'Iterative') %>%
#   rename(b = Estimate, se = 'Std.Error', t.stat = 't value', p = 'Pr(>|t|)')
# radial.ivw.exact_f <- as_tibble(radial.out$fe.coef) %>%
#   mutate(model = 'Exact (FE)') %>%
#   rename(b = Estimate, se = 'Std.Error', t.stat = 't value', p = 'Pr(>|t|)')
# radial.ivw.exact_r <- as_tibble(radial.out$re.coef) %>%
#   mutate(model = 'Iterative') %>%
#   rename(b = Estimate, se = 'Std.Error', t.stat = 't value', p = 'Pr(>|t|)')

# radial.ivw %>%
#   bind_rows(radial.ivw.iterative) %>%
#   bind_rows(radial.ivw.exact_f) %>%
#   select(model, b, se, t.stat, p)

## Modified Q Statistic Output
heterogenity.out <- tibble(
  id.exposure = as.character(mrdat.out[1,'id.exposure']),
  id.outcome = as.character(mrdat.out[1,'id.outcome']),
  outcome = as.character(mrdat.out[1,'outcome']),
  exposure = as.character(mrdat.out[1,'exposure']),
  pt = mrdat.raw %>% slice(1) %>% pull(pt),
  outliers_removed = FALSE,
  n_outliers = sum(radial.out$data$Outliers == 'Outlier'),
  Q.statistic = radial.out$qstatistic,
  df = radial.out$df) %>%
  mutate(p = pchisq(Q.statistic, df, lower.tail = FALSE))

## If Outliers detected, calculate Modified Q statistic after outlier removal
if (sum(mrdat.out$Outliers, na.rm = T) > 0){
	message('\nRE-RUNNING RADIAL MR AFTER OUTLIERS REMOVED\n')
  radial_wo_outliers.dat <- (
	  mrdat.out %>%
	    filter(pval.outcome > 5e-8) %>%
	    filter(mr_keep == TRUE) %>%
	    filter(Outliers == FALSE) %>%
	    dat_to_RadialMR(.)
	)[[1]]

  sink(paste0(out, '.txt'), append=TRUE, split=FALSE)
  radial_wo_outliers.out <- ivw_radial(radial_wo_outliers.dat, weights = 3, alpha = 0.05/nrow(radial_wo_outliers.dat))
  sink()

  heterogenity_wo_outliers.out <- tibble(
    id.exposure = as.character(mrdat.out[1,'id.exposure']),
    id.outcome = as.character(mrdat.out[1,'id.outcome']),
    outcome = as.character(mrdat.out[1,'outcome']),
    exposure = as.character(mrdat.out[1,'exposure']),
    pt = mrdat.raw %>% slice(1) %>% pull(pt),
    outliers_removed = TRUE,
    n_outliers = sum(radial_wo_outliers.out$data$Outliers == 'Outlier'),
    Q.statistic = radial_wo_outliers.out$qstatistic,
    df = radial_wo_outliers.out$df) %>%
    mutate(p = pchisq(Q.statistic, df, lower.tail = FALSE))
  heterogenity.out <- heterogenity.out %>% bind_rows(heterogenity_wo_outliers.out)
}

message('\nWRITING OUT RESULTS\n')
write_tsv(heterogenity.out, paste0(out, '_modifiedQ.txt'))
write_csv(mrdat.out, paste0(out, '_MRdatRadial.csv'))
