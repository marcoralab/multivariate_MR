#!/usr/bin/Rscript
message("Begining Harmonization \n")
### ===== Command Line Arguments ===== ##
exposure.summary = snakemake@input[["ExposureSummary"]] # Exposure summary statistics
outcome.summary = snakemake@input[["OutcomeSummary"]] # Outcome Summary statistics
proxy.snps = snakemake@input[["ProxySNPs"]]
pt = snakemake@params[["Pthreshold"]]
ExposureCode = snakemake@params[["excposurecode"]]
OutcomeCode = snakemake@params[["outcomecode"]]
out.harmonized = snakemake@output[["Harmonized"]] # SPECIFY THE OUTPUT FILE

### ===== Load packages ===== ###
#suppressMessages(library(plyr))
library(dplyr)  ## For data wrangling
library(readr)
library(TwoSampleMR) ## For conducting MR https://mrcieu.github.io/TwoSampleMR/

### ===== Read In Data ===== ###
message("READING IN EXPOSURE \n")
exposure.dat <- read_tsv(exposure.summary)

message("READING IN OUTCOME \n")
outcome.dat <- read_tsv(outcome.summary)

message("READING IN PROXY SNPs \n")
proxy.dat <- read_csv(proxy.snps) %>%
  filter(proxy.outcome == TRUE) %>%
  select(proxy.outcome, target_snp, proxy_snp, ALT, REF, ALT.proxy, REF.proxy) %>%
  mutate(SNP = target_snp) %>%
  rename(target_snp.outcome = target_snp, proxy_snp.outcome = proxy_snp, target_a1.outcome = ALT, target_a2.outcome = REF, proxy_a1.outcome = ALT.proxy, proxy_a2.outcome = REF.proxy)


### ===== Harmonization ===== ###
message("Harmonizing Exposure and Outcome \n")
mr_exposure.dat <- format_data(exposure.dat, type = 'exposure',
                            snp_col = 'SNP',
                            beta_col = "BETA",
                            se_col = "SE",
                            eaf_col = "AF",
                            effect_allele_col = "ALT",
                            other_allele_col = "REF",
                            pval_col = "P",
                            z_col = "Z",
                            samplesize_col = "N",
                            ncase_col = "N_CASES",
                        ncontrol_col = "N_CTRLS",
                            phenotype_col = 'TRAIT',
                            chr_col = 'CHROM',
                            pos_col = 'POS')
mr_exposure.dat$exposure =  ExposureCode

# Format LOAD
mr_outcome.dat <- format_data(outcome.dat, type = 'outcome',
                                 snp_col = 'SNP',
                                 beta_col = "BETA",
                                 se_col = "SE",
                                 eaf_col = "AF",
                                 effect_allele_col = "ALT",
                                 other_allele_col = "REF",
                                 pval_col = "P",
                                z_col = "Z",
                                samplesize_col = "N",
                            	ncase_col = "N_CASES",
                        	ncontrol_col = "N_CTRLS",
                                phenotype_col = 'TRAIT',
                                chr_col = 'CHROM',
                            	pos_col = 'POS')
mr_outcome.dat$outcome =  OutcomeCode

if(plyr::empty(proxy.dat) == FALSE){
  mr_outcome.dat <- left_join(mr_outcome.dat, proxy.dat, by = 'SNP')
}


# harmonize LOAD
harmonized.MRdat <- harmonise_data(mr_exposure.dat, mr_outcome.dat) %>%
  as_tibble() %>%
  mutate(pt = pt)

## Write out Harmonized data
write_csv(harmonized.MRdat, out.harmonized)
