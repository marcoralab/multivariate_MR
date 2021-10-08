if(any(grepl("conda", .libPaths(), fixed = TRUE))){
  message("Setting libPaths")
  df = .libPaths()
  conda_i = which(grepl("conda", df, fixed = TRUE))
  .libPaths(c(df[conda_i], df[-conda_i]))
}
.libPaths(c(.libPaths(), "/hpc/users/harern01/.Rlib"))

message("load in package \n")
library(tidyr)
library(readr)
library(dplyr)
library(MVMR)

message(" Command Line Arguments \n")
### ===== Command Line Arguments ===== ##
bmi.ss = snakemake@input[["bmi"]] # Exposure summary statistics
bmi.SNP = snakemake@params[["bmi_snplists"]] # significant BMI SNPs
t2d.ss = snakemake@input[["t2d"]] # Exposure Summary statistics
t2d.SNP = snakemake@params[["t2d_snplists"]] # significant T2D SNPs
COVID.ss = snakemake@input[["covid"]] # Outcome Summary statistics
outfile1 = snakemake@output[["outfile1"]] # Output
outfile2 = snakemake@output[["outfile2"]] # Output
outfile3 = snakemake@output[["outfile3"]] # Output
outfile4 = snakemake@output[["outfile4"]] # Output

# ========================================================================================== ##
##                                 Format exposure SNPlist

bmi.snps <- read_table2(bmi.SNP) %>%
  select(SNP,BETA,SE,TRAIT, P )%>%
  rename(BMI_beta = BETA, BMI_se = SE)

t2d.snps <- read_table2(t2d.SNP)%>%
  select(SNP,BETA,SE,TRAIT, P )%>%
  rename(T2D_beta = BETA, T2D_se = SE)

snps <- bind_rows(bmi.snps, t2d.snps) %>%
  as_tibble() %>%
  distinct(SNP) %>%
  pull(SNP)

# ========================================================================================== ##
##                                 Load exposure summary statistics data

bmi <- read_tsv(bmi.ss) %>%
    select(SNP,BETA,SE,TRAIT,P)
t2d <- read_tsv(t2d.ss)%>%
    select(SNP,BETA,SE,TRAIT,P)

# ========================================================================================== ##
##                                 Select SNPs from each exposure trait

bmi_t2d_snps <- bmi %>% filter(SNP %in% snps)%>%
    rename(BMI_beta = BETA, BMI_se = SE)

t2d_bmi_snps <- t2d %>% filter(SNP %in% snps)%>%
      rename(T2D_beta = BETA, T2D_se = SE)

## Full join all the exposure
exposure_dat <-full_join(bmi_t2d_snps,t2d_bmi_snps, by = "SNP" )%>%
        select(SNP,BMI_beta,BMI_se,T2D_beta, T2D_se )


##                                      COVID-19
# ========================================================================================== ##
##                                      OUTCOME: COVID
##  Read in covid outcome Files
covid <- read_tsv(COVID.ss)%>%
          select(SNP,BETA,SE,TRAIT,P)

covid_bmi_t2d_snps <- covid %>% filter(SNP %in% snps)%>%
  rename(COVID_beta = BETA, COVID_se = SE)

covid_mvmr <- full_join(exposure_dat,covid_bmi_t2d_snps, by = "SNP" )%>%
  select(BMI_beta,BMI_se,T2D_beta, T2D_se, COVID_beta ,COVID_se, SNP ) %>%
  na.omit(covid_mvmr)

## Format summary data
covid_mvmr.data <- format_mvmr(BXGs = covid_mvmr[,c(1,3)],
                        BYG = covid_mvmr[,5],
                        seBXGs = covid_mvmr[,c(2,4)],
                        seBYG = covid_mvmr[,6],
                        RSID = covid_mvmr[,7])

## Test for weak instruments
covid_mvmrcovmatrix <- matrix(c(1,0.239,0.239,1), nrow = 2, ncol = 2)
covid_Xcovmat <- phenocov_mvmr(covid_mvmrcovmatrix,covid_mvmr.data[,6:7])
covid_stre <- strength_mvmr(r_input = covid_mvmr.data, gencov = covid_Xcovmat)
covid_stre <- data.frame(covid_stre)

## Test for horizontal pleiotropy using conventional Q-statistic estimation

covid_pres <- pleiotropy_mvmr(r_input = covid_mvmr.data, gencov = covid_Xcovmat)
covid_pres <- data.frame(covid_pres$Qstat, covid_pres$Qpva)
covid_pres  <- covid_pres %>% rename(Qstat= covid_pres.Qstat, Qpva = covid_pres.Qpva)

## Estimate causal effects
res_covid <- ivw_mvmr(r_input = covid_mvmr.data)
exposure1_lo_ci <- res_covid[1] - 1.96 *res_covid[3]
exposure2_lo_ci <- res_covid[2] - 1.96 *res_covid[4]
exposure1_up_ci <- res_covid[1] + 1.96 *res_covid[3]
exposure2_up_ci <- res_covid[2] + 1.96 *res_covid[4]
exposure1_or = exp(res_covid[1])
exposure2_or = exp(res_covid[2])
exposure1_or_lci95 = exp(exposure1_lo_ci)
exposure2_or_lci95 = exp(exposure2_lo_ci)
exposure1_or_uci95 = exp(exposure1_up_ci)
exposure2_or_uci95 = exp(exposure2_up_ci)

res_covid_OR <- data.frame(exposure1_or,exposure2_or,exposure1_or_lci95,exposure2_or_lci95,exposure1_or_uci95,exposure2_or_uci95)%>%
  rename(BMI_OR = exposure1_or, T2D_OR = exposure2_or, BMI_OR_lci95 = exposure1_or_lci95, T2D_OR_lci95 = exposure2_or_lci95, BMI_OR_uci95 = exposure1_or_uci95, T2D_OR_uci95 = exposure2_or_uci95)

## Robust causal effect estimation

robust_res_covid<- qhet_mvmr(covid_mvmr.data, covid_mvmrcovmatrix, CI = F, iterations = 100)
robust_res_covid <- data.frame(robust_res_covid)

## ================= Write Out ================= ##

# instrument strength
covid_stre %>% write_tsv(outfile1)

# pleiotropy out
covid_pres %>% write_tsv(outfile2)

# Results out
res_covid_OR %>% write_tsv(outfile3)

## Robust causal effect estimation out
robust_res_covid %>% write_tsv(outfile4)
