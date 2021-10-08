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

`%nin%` = Negate(`%in%`)

## =================================== ##
# setwd('~/GitCode/mvmr_hight_ICV_AD')
# bmi.name = "Jiang2019height"
# bmi.name = "NealeLab2018height"
# bmi.name = "CanelaXandri2018height"
# diabetes.name = "Andrews2019hcicvbv"
# diabetes.name = "Haworth2019hcicv"
# COVID_B2__EUR.name = "Kunkle2019load_stage123"
# out =

### ===== Command Line Arguments ===== ##
args = commandArgs(trailingOnly = TRUE) # Set arguments from the command line
bmi.name = snakemake@params[["mvexp1"]] # Exposure summary statistics
diabetes.name = snakemake@params[["mvexp2"]] # Outcome Summary statistics
COVID_B2__EUR.name = snakemake@params[["mvout"]]
snplist.path = snakemake@input[["mvmr_snplists"]]
out.path = snakemake@output[["MVMRresults"]]

# Load Height GWAS summary stats
bmi.raw <- read_tsv(paste0('../MRcovid/data/formated/Yengo2018bmi/', bmi.name, '_formated.txt.gz'), comment = '##',
                       col_types = list(SNP = col_character(),
                                        CHROM = col_double(),
                                        POS = col_double(),
                                        REF = col_character(),
                                        AF = col_double(),
                                        TRAIT = col_character(),
                                        BETA= col_double(),
                                        SE = col_double(),
                                        Z = col_double(),
                                        P = col_double(),
                                        N = col_double())) %>%
  mutate(TRAIT = bmi.name)%>%
  filter(AF != 0)

# Load Brain volume Summary Stats
diabetes.raw <- read_tsv(paste0('../MRcovid/data/formated/Mahajan2018t2d/', diabetes.name, '_formated.txt.gz'), comment = '##',
                       col_types = list(SNP = col_character(),
                                        CHROM = col_double(),
                                        POS = col_double(),
                                        REF = col_character(),
                                        AF = col_double(),
                                        TRAIT = col_character(),
                                        BETA= col_double(),
                                        SE = col_double(),
                                        Z = col_double(),
                                        P = col_double(),
                                        N = col_double())) %>%
  mutate(TRAIT = diabetes.name) %>%
  filter(AF != 0)

# Load Brain volume Summary Stats
COVID_B2__EUR.raw <- read_tsv(paste0('../MRcovid/data/formated/covidhgi2020B2v6alleur/', COVID_B2__EUR.name, '_formated.txt.gz'), comment = '##',
                     col_types = list(SNP = col_character(),
                                      CHROM = col_double(),
                                      POS = col_double(),
                                      REF = col_character(),
                                      AF = col_double(),
                                      TRAIT = col_character(),
                                      BETA= col_double(),
                                      SE = col_double(),
                                      Z = col_double(),
                                      P = col_double(),
                                      N = col_double())) %>%
  mutate(TRAIT = COVID_B2__EUR.name) %>%
  format_data(., type = 'outcome',
              snp_col = 'SNP',
              beta_col = "BETA",
              se_col = "SE",
              eaf_col = "AF",
              effect_allele_col = "ALT",
              other_allele_col = "REF",
              pval_col = "P",
              z_col = "Z",
              samplesize_col = "N",
              phenotype_col = 'TRAIT',
              chr_col = 'CHROM',
              pos_col = 'POS'
  ) %>%
  as_tibble()


# ## ========================================================================================== ##
# ##                                      Height -> Load
# height_load.path = paste0('input/', bmi.name, '/', COVID_B2__EUR.name, '/', bmi.name, '_5e-8_', COVID_B2__EUR.name, '_MRdatRadial.csv')
# height_load.dat <- read_csv(height_load.path)
#
# height_ad.res <- height_load.dat %>%
#   filter(mr_keep == TRUE, Outliers == FALSE) %>%
#   mr(., method_list = c("mr_ivw_fe", "mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
# generate_odds_ratios(height_ad.res)
#
# ## ========================================================================================== ##
# ##                                      Height -> Load
# height_bv.path = paste0('input/', bmi.name, '/', diabetes.name, '/', bmi.name, '_5e-8_', diabetes.name, '_MRdatRadial.csv')
# height_bv.dat <- read_csv(height_bv.path)
#
# height_hcicv.res <- height_bv.dat  %>%
#   filter(mr_keep == TRUE, Outliers == FALSE) %>%
#   mr(., method_list = c("mr_ivw_fe", "mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
# height_hcicv.res
#
# ## ========================================================================================== ##
# ##                                      BV -> Load
# bv_load.path = paste0('input/', diabetes.name, '/', COVID_B2__EUR.name, '/', diabetes.name, '_5e-8_', COVID_B2__EUR.name, '_MRdatRadial.csv')
# bv_load.dat <- read_csv(bv_load.path)
#
# bv_load.res <- bv_load.dat %>%
#   filter(mr_keep == TRUE, Outliers == FALSE) %>%
#   mr(., method_list = c("mr_ivw_fe", "mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
# generate_odds_ratios(bv_load.res)
#
# ## ========================================================================================== ##
# ##                                      bv -> Height
# bv_height.path = paste0('input/', diabetes.name, '/', bmi.name, '/', diabetes.name, '_5e-8_', bmi.name, '_MRdatRadial.csv')
# bv_height.dat <- read_csv(bv_height.path)
#
# bv_height.res <- bv_height.dat %>%
#   filter(mr_keep == TRUE, Outliers == FALSE) %>%
#   mr(., method_list = c("mr_ivw_fe", "mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
# bv_height.res
#
# ## ========================================================================================== ##
# ##                                      MVMR

snps <- read_table2(snplist.path)

## Select height + HC+ICV SNPs from HC+ICV GWAS
## Clump result dataframe
diabetes_all <- diabetes.raw %>%
  filter(SNP %in% snps$SNP) %>%
  format_data(., type = 'exposure',
              snp_col = 'SNP',
              beta_col = "BETA",
              se_col = "SE",
              eaf_col = "AF",
              effect_allele_col = "ALT",
              other_allele_col = "REF",
              pval_col = "P",
              z_col = "Z",
              samplesize_col = "N",
              phenotype_col = 'TRAIT',
              chr_col = 'CHROM',
              pos_col = 'POS')  %>%
  as_tibble()

## Extract clumped SNPs from Height GWAS
bmi_all <- bmi.raw %>%
  filter(SNP %in% snps$SNP) %>%
  format_data(., type = 'exposure',
              snp_col = 'SNP',
              beta_col = "BETA",
              se_col = "SE",
              eaf_col = "AF",
              effect_allele_col = "ALT",
              other_allele_col = "REF",
              pval_col = "P",
              z_col = "Z",
              samplesize_col = "N",
              phenotype_col = 'TRAIT',
              chr_col = 'CHROM',
              pos_col = 'POS') %>%
  as_tibble()

## Bind GWAS togther
exposure_dat <- bind_rows(diabetes_all, bmi_all)

## harmonize exposures and outcomes
mvdat <- mv_harmonise_data(exposure_dat, COVID_B2__EUR.raw, harmonise_strictness = 1)
mvmr_ad <- mv_multiple(mvdat, pval_threshold = 5e-08)

write_tsv(mvmr_ad$result, out.path)

#
# harmonise_data(bmi_all, COVID_B2__EUR.raw) %>%
#   filter(pval.exposure < 5e-8) %>%
#   mr(., method_list = c("mr_ivw_fe", "mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
#
# harmonise_data(diabetes_all, COVID_B2__EUR.raw) %>%
#   filter(pval.exposure < 5e-8) %>%
#   mr(., method_list = c("mr_ivw_fe", "mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
#
