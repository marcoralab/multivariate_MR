library(readxl)
library(tidyverse)
library(ggman)
library(TwoSampleMR)
z2se <- function(p, n, z){1/sqrt(2*p*(1-p)*(n+z^2))}
setwd('/Users/sheaandrews/GitCode/mvmr_hight_ICV_AD/')

## Study Overlap
kunkle <- read_excel('/Users/sheaandrews/GitCode/mvmr_hight_ICV_AD/data/StudyOverLap.xlsx', sheet = 'Kunkle') %>% 
  mutate(study = 'kunkle', trait = 'AD') %>% 
  setNames(paste0('kunkle.', names(.)))

adams <- read_excel('/Users/sheaandrews/GitCode/mvmr_hight_ICV_AD/data/StudyOverLap.xlsx', sheet = 'Adams') %>% 
  mutate(study = 'adams', trait = 'ICV')
haworth <- read_excel('/Users/sheaandrews/GitCode/mvmr_hight_ICV_AD/data/StudyOverLap.xlsx', sheet = 'Haworth') %>% 
  mutate(study = 'haworth', trait = 'HC')

andrews = bind_rows(adams, haworth) %>% 
  add_row(Consortium = 'UKBB', cohort = 'UKBB', n = 8428, study = 'Elliott', trait = 'BV') %>% 
  setNames(paste0('andrews.', names(.)))
wood <- read_excel('/Users/sheaandrews/GitCode/mvmr_hight_ICV_AD/data/StudyOverLap.xlsx', sheet = 'Wood') %>% 
  mutate(study = 'wood', trait = 'Height') %>% 
  setNames(paste0('wood.', names(.)))
watanabe <- tibble(Consortium = 'UKBB', cohort = 'UKBB', n = 385748, study = 'Watanabe', trait = 'Height') %>% 
  setNames(paste0('watanabe.', names(.)))


study = full_join(kunkle, andrews, by = c('kunkle.cohort' = 'andrews.cohort')) %>% 
  full_join(watanabe, by = c('kunkle.cohort' = 'watanabe.cohort')) %>% 
  rename(cohort = kunkle.cohort) %>% 
  select(cohort, kunkle.n, andrews.n, watanabe.n) %>% 
  mutate(overlap_kunkle_andrews = ifelse(kunkle.n < andrews.n, kunkle.n, andrews.n), 
         overlap_watanabe_andrews = ifelse(watanabe.n < andrews.n, watanabe.n, andrews.n))
print(study, n = Inf)

summarise(study, 
          andrews.n = sum(andrews.n, na.rm = T),
          kunkle.n = sum(kunkle.n, na.rm = T), 
          adams.n = sum(andrews.n, na.rm = T), 
          watanabe.n = sum(watanabe.n, na.rm = T),
          overlap_kunkle_andrews = sum(overlap_kunkle_andrews, na.rm = T), 
          pc_kunkle_andrews = (overlap_kunkle_andrews / kunkle.n) * 100, 
          pc_andrews_kunkle = (overlap_kunkle_andrews / andrews.n) * 100, 
          overlap_watanabe_andrews = sum(overlap_watanabe_andrews, na.rm = T), 
          pc_watanabe_andrews = (overlap_watanabe_andrews / watanabe.n) * 100)

## Brain Volume
bv.raw <- read_tsv('data/GWAS/Elliott2018bv009.txt.gz')

count(bv, P < 5e-8)
ggman(filter(bv, P < 0.05), snp = 'RSID', bp = 'POS', chrom = 'CHROM', pvalue = 'P', ymax = 15)

bv.dat <- bv %>% 
  filter(P < 5e-8) %>% 
  format_data(., type = 'exposure',
              snp_col = 'RSID',
              beta_col = "BETA",
              se_col = "SEBETA",
              eaf_col = "EAF",
              effect_allele_col = "ALT",
              other_allele_col = "REF",
              pval_col = "P",
              z_col = "Z",
              samplesize_col = "N",
              phenotype_col = 'TRAIT') %>%
  clump_data(.) %>% 
  as_tibble() 

## HC 

hc.raw <- read_table2('data/GWAS/Haworth2019hc.chrall.CPRA_b37.tsv.gz', comment = '##')
hc.raw %>% 
  filter(P < 5e-8) %>% 
  mutate(TRAIT = 'icvhcbc') %>% 
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

## ICV 

icv.raw <- read_tsv('data/GWAS/Adams2014icv_fixed.gz')

## ICV + HC

hcicv.raw <- read_table2('data/GWAS/Haworth2019hcicv.chrall.CPRA_b37.tsv.gz', comment = '##')

hcicv.raw %>% 
  filter(P < 5e-8) %>% 
  mutate(TRAIT = 'icvhcbc') %>% 
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

## HC + ICV + BC
hcicvbc.raw <- read_tsv('/Users/sheaandrews/LOAD_minerva/dummy/shea/Projects/mvmr_hight_ICV_AD/data/GWAS/hcicvbc.metal1.tbl')

out <- hcicvbc.raw %>%
  rename(P = `P-value`, EA = Allele1, OA = Allele2, EAF = Freq1, SNP = MarkerName, N = Weight, CHROM = Chromosome, POS = Position) %>%
  mutate(EA = toupper(EA), 
         OA = toupper(OA)) %>%
  mutate(SE = z2se(EAF, N, Zscore)) %>%
  mutate(Beta = Zscore*SE) %>% 
  arrange(CHROM, POS)

write_tsv(out, gzfile('data/GWAS/andrews2019hcicvbc.txt.gz'))

out %>% 
  filter(!is.na(CHROM)) %>% 
  filter(P < 0.05) %>% 
  ggman(., snp = 'SNP', chrom = 'CHROM', bp = 'POS', pvalue = 'P', ymin = 0, ymax = 25, relative.positions = TRUE) + theme_bw()

icvhcbc.dat <- out %>% 
  filter(P < 5e-8) %>% 
  mutate(TRAIT = 'icvhcbc') %>% 
  format_data(., type = 'exposure',
              snp_col = 'SNP',
              beta_col = "Beta",
              se_col = "SE",
              eaf_col = "EAF",
              effect_allele_col = "EA",
              other_allele_col = "OA",
              pval_col = "P",
              z_col = "Zscore",
              samplesize_col = "N",
              phenotype_col = 'TRAIT', 
              chr_col = 'CHROM', 
              pos_col = 'POS') %>%
  clump_data(.) %>% 
  as_tibble() 

## Height - Wood 


count(wood, p < 5e-8)
ggman(filter(wood, p < 0.05), snp = 'MarkerName', bp = 'POS', chrom = 'CHROM', pvalue = 'p', ymax = 30)

height.dat <- wood %>% 
  filter(p < 5e-8) %>% 
  format_data(., type = 'exposure',
              snp_col = 'MarkerName',
              beta_col = "b",
              se_col = "SE",
              eaf_col = "Freq.Allele1.HapMapCEU",
              effect_allele_col = "Allele1",
              other_allele_col = "Allele2",
              pval_col = "p",
              z_col = "Z",
              samplesize_col = "N",
              phenotype_col = 'TRAIT') %>%
  clump_data(.) %>% 
  as_tibble() 






