library(tidyverse)
library(TwoSampleMR)
library(RadialMR)
`%nin%` = Negate(`%in%`)

std_beta = function(z, eaf, n){
  std.b = z/sqrt(2 * eaf * (1 - eaf) * (n + z^2))
  std.b
}

std_se = function(z, eaf, n){
  std.se = 1/sqrt(2 * eaf * (1 - eaf) * (n + z^2))
  std.se
}

## ========================================================================================== ##
# setwd('~/GitCode/mvmr_hight_ICV_AD')
# height.name = "Jiang2019height"
# bv.name = "Andrews2019hcicvbv"
# load.name = "Kunkle2019load_stage123"
# out =

### ===== Command Line Arguments ===== ##
args = commandArgs(trailingOnly = TRUE) # Set arguments from the command line
height.name = args[1] # Exposure summary statistics
bv.name = args[2] # Outcome Summary statistics
load.name = args[3]
out = args[3]


# Load Height GWAS summary stats
NealLab.raw <- read_tsv('data/GWAS/NealeLab2018height.chrall.CPRA_b37.tsv.gz', comment = '##') %>%
  mutate(TRAIT = 'Height_NealLab') %>%
  select(-N_CASES, -N_CTRLS, -OR, -OR_L95, -OR_U95, -DIR, -P_HET, -G1000_ID, -G1000_VARIANT, -DBSNP_VARIANT, -OLD_ID, -OLD_VARIANT)
watanabe.raw <- read_table2('data/GWAS/Watanabe2019height.chrall.CPRA_b37.tsv.gz', comment = '##') %>%
  mutate(TRAIT = 'Height_Watanabe') %>%
  select(-N_CASES, -N_CTRLS, -OR, -OR_L95, -OR_U95, -DIR, -P_HET, -G1000_ID, -G1000_VARIANT, -DBSNP_VARIANT, -OLD_ID, -OLD_VARIANT)
CanelaXandri.raw <- read_table2('data/GWAS/CanelaXandri2018height.chrall.CPRA_b37.tsv.gz', comment = '##') %>%
  mutate(TRAIT = 'Height_CanelaXandri') %>%
  select(-N_CASES, -N_CTRLS, -OR, -OR_L95, -OR_U95, -DIR, -P_HET, -G1000_ID, -G1000_VARIANT, -DBSNP_VARIANT, -OLD_ID, -OLD_VARIANT)
Jiang.raw <- read_table2('data/GWAS/Jiang2019height.chrall.CPRA_b37.tsv.gz', comment = '##') %>%
  mutate(TRAIT = 'Height_Jiang') %>%
  select(-N_CASES, -N_CTRLS, -OR, -OR_L95, -OR_U95, -DIR, -P_HET, -G1000_ID, -G1000_VARIANT, -DBSNP_VARIANT, -OLD_ID, -OLD_VARIANT)
wood.raw <- read_tsv('data/GWAS/Wood2014height.chrall.CPRA_b37.tsv.gz', comment = '##') %>%
  mutate(TRAIT = 'Height_wood') %>%
  select(-N_CASES, -N_CTRLS, -OR, -OR_L95, -OR_U95, -DIR, -P_HET, -G1000_ID, -G1000_VARIANT, -DBSNP_VARIANT, -OLD_ID, -OLD_VARIANT)
yengo.raw <- read_tsv('data/GWAS/Yengo2018height.chrall.CPRA_b37.tsv.gz', comment = '##')  %>%
  mutate(TRAIT = 'Height_Yengo') %>%
  select(-N_CASES, -N_CTRLS, -OR, -OR_L95, -OR_U95, -DIR, -P_HET, -G1000_ID, -G1000_VARIANT, -DBSNP_VARIANT, -OLD_ID, -OLD_VARIANT)


# Load Brain volume datasets
## Either use published HC + ICV OR meta-anlysis with BV
hcicv.raw <- read_tsv('data/GWAS/Haworth2019hcicv.chrall.CPRA_b37.tsv.gz', comment = '##')
hcicvdiabetes.raw <- read_tsv('data/GWAS/Andrews2019hcicvbv.chrall.CPRA_b37.tsv.gz', comment = '##') %>%
  filter(AF != 0)

# Load Kunkle AD GWAS Summary stats
fh.raw <- read_tsv('/Users/sheaandrews/Dropbox/Research/Data/Summary_Statisitics/FultonHoward2019load.chrall.CPRA_b37.tsv.gz', comment = '##')
kunkle.raw <- read_tsv('data/GWAS/Kunkle2019load_stage123.chrall.CPRA_b37.tsv.gz', comment = '##')

load.raw <- kunkle.raw %>%
  format_data(., type = 'outcome',
              snp_col = 'DBSNP_ID',
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


## ========================================================================================== ##
##                                      Height
height <- height.raw %>%
  filter(ALT %in% c('A', 'C', 'G', 'T')) %>%
  filter(REF %in% c('A', 'C', 'G', 'T')) %>%
  filter(P < 5e-8) %>%
  #filter(P < 5e-70) %>%
  mutate(sb = std_beta(Z, AF, N)) %>%
  mutate(sse = std_se(Z, AF, N)) %>%
  format_data(., type = 'exposure',
              snp_col = 'DBSNP_ID',
              #beta_col = "BETA",
              #se_col = "SE",
              beta_col = "sb",
              se_col = "sse",
              eaf_col = "AF",
              effect_allele_col = "ALT",
              other_allele_col = "REF",
              pval_col = "P",
              z_col = "Z",
              samplesize_col = "N",
              phenotype_col = 'TRAIT',
              chr_col = 'CHROM',
              pos_col = 'POS') %>%
#  clump_data(.) %>%
  group_split(chr.exposure, keep = TRUE) %>%
  map(., clump_data) %>%
  bind_rows() %>%
  as_tibble()

## height -> Load
height_load.dat <- harmonise_data(height, load.raw)

## Radial MR
height_load.radial <- (
  height_load.dat %>%
  filter(pval.outcome > 5e-8) %>%
  filter(mr_keep == TRUE) %>%
  dat_to_RadialMR(.)
  )[[1]]
height_load.radial.res <- ivw_radial(height_load.radial, weights = 3, alpha = 0.05/nrow(height_load.radial))
# plot_radial(height_load.radial.res)

## MR
height_load.mrdat <- height_load.dat %>%
  left_join(height_load.radial.res$data, by = 'SNP')

height_ad.res <- height_load.mrdat %>%
  filter(mr_keep == TRUE, Outliers == 'Variant') %>%
  mr(., method_list = c("mr_ivw_fe", "mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
generate_odds_ratios(height_ad.res)

## height -> HC + ICV
height_hcicv <- diabetes.raw %>%
  filter(DBSNP_ID %in% height$SNP) %>%
  format_data(., type = 'outcome',
              snp_col = 'DBSNP_ID',
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
              pos_col = 'POS')

height_hcicv.dat <- harmonise_data(height, height_hcicv)

## Radial MR
height_hcicv.radial <- (
  height_hcicv.dat %>%
    filter(pval.outcome > 5e-8) %>%
    filter(mr_keep == TRUE) %>%
    dat_to_RadialMR(.)
)[[1]]

height_hcicv.radial.res <- ivw_radial(height_hcicv.radial, weights = 3, alpha = 0.05/nrow(height_hcicv.radial))
# plot_radial(height_hcicv.radial.res)
## MR
height_hcicv.mrdat <- height_hcicv.dat %>%
  left_join(height_hcicv.radial.res$data, by = 'SNP')

height_hcicv.res <- height_hcicv.mrdat %>%
  filter(mr_keep == TRUE, Outliers == 'Variant') %>%
  mr(., method_list = c("mr_ivw_fe", "mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
height_hcicv.res

## ========================================================================================== ##
##                                      HC + ICV
bv <- diabetes.raw %>%
#  filter(P < 5e-6) %>%
  filter(P < 5e-8) %>%
  format_data(., type = 'exposure',
              snp_col = 'DBSNP_ID',
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
  clump_data(.) %>%
  as_tibble(.)

# BV -> Load
## harmonize
bv_load.dat <- harmonise_data(bv, load.raw)

## Radial MR
bv_load.radial <- (
  bv_load.dat %>%
    filter(pval.outcome > 5e-8) %>%
    filter(mr_keep == TRUE) %>%
    dat_to_RadialMR(.)
)[[1]]

bv_load.radial.res <- ivw_radial(bv_load.radial, weights = 3, alpha = 0.05/nrow(bv_load.radial))
plot_radial(bv_load.radial.res)
## MR
bv_load.mrdat <- bv_load.dat %>%
  left_join(bv_load.radial.res$data, by = 'SNP')

bv_load.res <- bv_load.mrdat %>%
  filter(mr_keep == TRUE, Outliers == 'Variant') %>%
  mr(., method_list = c("mr_ivw_fe", "mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
generate_odds_ratios(bv_load.res)


## bv -> Height
bv_height.snps <- height.raw %>%
  filter(DBSNP_ID %in% bv$SNP) %>%
  mutate(TRAIT = 'Height') %>%
  format_data(., type = 'outcome',
              snp_col = 'DBSNP_ID',
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
              pos_col = 'POS')

bv_height.dat <- harmonise_data(bv, bv_height.snps)

## Radial MR
bv_height.radial <- (
  bv_height.dat %>%
    filter(pval.outcome > 5e-8) %>%
    filter(mr_keep == TRUE) %>%
    dat_to_RadialMR(.)
)[[1]]

bv_height.radial.res <- ivw_radial(bv_height.radial, weights = 3, alpha = 0.05/nrow(bv_height.radial))
plot_radial(bv_height.radial.res, radial_scale = F)

## MR
bv_height.mrdat <- bv_height.dat %>%
  left_join(bv_height.radial.res$data, by = 'SNP')

bv_height.res <- bv_height.mrdat %>%
  filter(mr_keep == TRUE, Outliers == 'Variant') %>%
  mr(., method_list = c("mr_ivw_fe", "mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
bv_height.res

bv_height.dat %>%
  left_join(bv_height.radial.res$data, by = 'SNP') %>%
  mr_heterogeneity(.)
directionality_test(bv_height.mrdat)

## ========================================================================================== ##
##                                      MVMR

bv.outliers <- bv_load.mrdat %>%
  filter(Outliers == 'Outlier') %>%
  mutate(SNP = as.character(SNP)) %>%
  pull(SNP)

height.outliers <- height_load.mrdat %>%
  filter(Outliers == 'Outlier') %>%
  mutate(SNP = as.character(SNP)) %>%
  pull(SNP)

filter(bv, SNP %in% height.outliers)
filter(height, SNP %in% bv.outliers)

snps <- bind_rows(height, bv) %>%
  as_tibble() %>%
  distinct(SNP) %>%
  filter(SNP %nin% bv.outliers) %>%
  filter(SNP %nin% height.outliers) %>%
  pull(SNP)

## Select height + HC+ICV SNPs from HC+ICV GWAS
## Clump result dataframe
diabetes_all <- diabetes.raw %>%
  filter(DBSNP_ID %in% snps) %>%
  format_data(., type = 'exposure',
              snp_col = 'DBSNP_ID',
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
  clump_data(.) %>%
  as_tibble()

## Extract clumped SNPs from Height GWAS
height_all <- height.raw %>%
  filter(DBSNP_ID %in% diabetes_all$SNP) %>%
  format_data(., type = 'exposure',
              snp_col = 'DBSNP_ID',
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
exposure_dat <- bind_rows(diabetes_all, height_all)

## harmonize exposures and outcomes
mvdat <- mv_harmonise_data(exposure_dat, load.raw, harmonise_strictness = 1)
mvmr_ad <- mv_multiple(mvdat, pval_threshold = 5e-08)
mvmr_ad

id_exposure <- c(299, 300, 302)
id_outcome <- 7

exposure_dat <- mv_extract_exposures(id_exposure)

bv_load <- harmonise_data(bv, load.raw)
height_load <- harmonise_data(height, load.raw)
mr(bv_load)
mr(height_load)
