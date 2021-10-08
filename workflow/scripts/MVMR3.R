setwd('/sc/arion/projects/LOAD/harern01/projects/')
library(tidyverse)
library(glue)
library(TwoSampleMR)
library(remotes)
library (MVMR)
library(RadialMR)
`%nin%` = Negate(`%in%`)

bmi_ss.path = glue("MRcovid/data/formated/Yengo2018bmi/Yengo2018bmi_formated.txt.gz")
diabetes_ss.path = glue("MRcovid/data/formated/Mahajan2018t2d/Mahajan2018t2d_formated.txt.gz")
Yengo2018bmi_5e_8_SNPs = glue("mvmr_bmi_T2D_COVID/input/Yengo2018bmi/Yengo2018bmi_5e-8_SNPs.txt")
Mahajan2018t2d_5e_8_SNP = glue("mvmr_bmi_T2D_COVID/input/Mahajan2018t2d/Mahajan2018t2d_5e-8_SNPs.txt")
covidA2_ss.path = glue("MRcovid/data/formated/covidhgi2020A2v6alleur/covidhgi2020A2v6alleur_formated.txt.gz")
covidB2_ss.path = glue("MRcovid/data/formated/covidhgi2020B2v6alleur/covidhgi2020B2v6alleur_formated.txt.gz")
covidC2_ss.path = glue("MRcovid/data/formated/covidhgi2020C2v6alleur/covidhgi2020C2v6alleur_formated.txt.gz")

# ========================================================================================== ## 
##                                      EXPOSURE

##  Read SNPs that are associated w/ exposures at p < 5e-8
Yengo2018bmi_5e_8_SNPs <- read_table2(Yengo2018bmi_5e_8_SNPs) %>% 
  select(SNP,BETA,SE,TRAIT, P )%>% 
  rename(BMI_beta = BETA, BMI_se = SE)

Mahajan2018t2d_5e_8_SNP <- read_table2(Mahajan2018t2d_5e_8_SNP)%>% 
  select(SNP,BETA,SE,TRAIT, P )%>% 
  rename(T2D_beta = BETA, T2D_se = SE)

snps <- bind_rows(Yengo2018bmi_5e_8_SNPs, Mahajan2018t2d_5e_8_SNP) %>% 
  as_tibble() %>% 
  distinct(SNP) %>% 
  pull(SNP) 

##  Read in BMI Exposure Files

bmi_ss <- read_tsv(bmi_ss.path)%>% 
  select(SNP,BETA,SE,TRAIT, P ) 

bmi_t2d_snps <- bmi_ss %>% filter(SNP %in% snps)%>% 
  rename(BMI_beta = BETA, BMI_se = SE)
  

##  Read in T2D Exposure Files
diabetes_ss <- read_tsv(diabetes_ss.path)%>% 
  select(SNP,BETA,SE,TRAIT, P ) 

t2d_bmi_snps <- diabetes_ss %>% filter(SNP %in% snps)%>% 
  rename(T2D_beta = BETA, T2D_se = SE)

## Full join all the exposures
exposure_mvmr <-full_join(bmi_t2d_snps,t2d_bmi_snps, by = "SNP" )%>%
  select(SNP,BMI_beta,BMI_se,T2D_beta, T2D_se ) 

##                                      COVID-A2
# ========================================================================================== ## 
##                                      OUTCOME
##  Read in covidA2 outcome Files
covidA2_ss <- read_tsv(covidA2_ss.path)%>% 
  select(SNP,BETA,SE,TRAIT, P ) 

covida2_bmi_t2d_snps <- covidA2_ss %>% filter(SNP %in% snps)%>% 
  rename(COVIDA2_beta = BETA, COVIDA2_se = SE)

covida2_mvmr <-full_join(exposure_mvmr,covida2_bmi_t2d_snps, by = "SNP" )%>%
  select(BMI_beta,BMI_se,T2D_beta, T2D_se, COVIDA2_beta ,COVIDA2_se, SNP ) 

covida2_mvmr <- na.omit(covida2_mvmr)


# ========================================================================================== ## 
##                                  Format summary data     

covida2_mvmr.data <- format_mvmr(BXGs = covida2_mvmr[,c(1,3)],
                      BYG = covida2_mvmr[,5],
                      seBXGs = covida2_mvmr[,c(2,4)],
                      seBYG = covida2_mvmr[,6],
                      RSID = covida2_mvmr[,7])

# ========================================================================================== ## 
##                                  Test for weak instruments 

sres <- strength_mvmr(r_input = covida2_mvmr.data, gencov = 0)

#Conditional F-statistics for instrument strength

             #exposure1  #exposure2
#F-statistic  30.63079   15.63525

mvmrcovmatrix<-matrix(c(1,-0.1,-0.05,-0.1,1,0.2,-0.05,0.2,1), nrow = 3, ncol = 3)
Xcovmat<-phenocov_mvmr(mvmrcovmatrix,covida2_mvmr.data[,5:7])
sres2 <- strength_mvmr(r_input = covida2_mvmr.data, gencov = Xcovmat)


# ========================================================================================== ## 
##            Test for horizontal pleiotropy using conventional Q-statistic estimation

pres <- pleiotropy_mvmr(r_input = covida2_mvmr.data, gencov = 0)

# ========================================================================================== ## 
##                              Estimate causal effects

res_covidA2 <- ivw_mvmr(r_input = covida2_mvmr.data)
exposure1_lo_ci <- res_covidA2[1] - 1.96 *res_covidA2[3]
exposure2_lo_ci <- res_covidA2[2] - 1.96 *res_covidA2[4]
exposure1_up_ci <- res_covidA2[1] + 1.96 *res_covidA2[3]
exposure2_up_ci <- res_covidA2[2] + 1.96 *res_covidA2[4]
exposure1_or = exp(res_covidA2[1])
exposure2_or = exp(res_covidA2[2])
exposure1_or_lci95 = exp(exposure1_lo_ci)
exposure2_or_lci95 = exp(exposure2_lo_ci)
exposure1_or_uci95 = exp(exposure1_up_ci)
exposure2_or_uci95 = exp(exposure2_up_ci)

res_covidA2_OR <- cbind(exposure1_or,exposure2_or,exposure1_or_lci95,exposure2_or_lci95,exposure1_or_uci95,exposure2_or_uci95)
res_covidA2_OR 

#Multivariable MR
#Estimate Std. Error  t value     Pr(>|t|)
#exposure1 0.31734282 0.07137197 4.446323 1.041169e-05
#exposure2 0.04059782 0.03379084 1.201445 2.300536e-01

# ========================================================================================== ## 
##                              Robust causal effect estimation

robust_res_covidA2<- qhet_mvmr(F.data, mvmrcovmatrix, CI = F, iterations = 100)

# ====================================================================================================================== ## 
# ====================================================================================================================== ## 
# ====================================================================================================================== ## 
# ====================================================================================================================== ##

##                                      COVID-B2
# ========================================================================================== ## 
##                                      OUTCOME
##  Read in covidB2 outcome Files
covidB2_ss <- read_tsv(covidB2_ss.path)%>% 
  select(SNP,BETA,SE,TRAIT, P ) 

covidb2_bmi_t2d_snps <- covidB2_ss %>% filter(SNP %in% snps)%>% 
  rename(COVIDB2_beta = BETA, COVIDB2_se = SE)

covidb2_mvmr <-full_join(exposure_mvmr,covidb2_bmi_t2d_snps, by = "SNP" )%>%
  select(BMI_beta,BMI_se,T2D_beta, T2D_se, COVIDB2_beta ,COVIDB2_se, SNP ) 

covidb2_mvmr <- na.omit(covidb2_mvmr)


# ========================================================================================== ## 
##                                  Format summary data     

covidb2_mvmr.data <- format_mvmr(BXGs = covidb2_mvmr[,c(1,3)],
                      BYG = covidb2_mvmr[,5],
                      seBXGs = covidb2_mvmr[,c(2,4)],
                      seBYG = covidb2_mvmr[,6],
                      RSID = covidb2_mvmr[,7])

# ========================================================================================== ## 
##                                  Test for weak instruments 

sres <- strength_mvmr(r_input = covidb2_mvmr.data, gencov = 0)

mvmrcovmatrix<-matrix(c(1,-0.1,-0.05,-0.1,1,0.2,-0.05,0.2,1), nrow = 3, ncol = 3)
Xcovmat<-phenocov_mvmr(mvmrcovmatrix,covidb2_mvmr.data[,5:7])
sres2 <- strength_mvmr(r_input = covidb2_mvmr.data, gencov = Xcovmat)


# ========================================================================================== ## 
##            Test for horizontal pleiotropy using conventional Q-statistic estimation

pres <- pleiotropy_mvmr(r_input = covidb2_mvmr.data, gencov = 0)
# Q-Statistic for instrument validity:


# ========================================================================================== ## 
##                              Estimate causal effects

res_covidB2 <- ivw_mvmr(r_input = covidb2_mvmr.data)
exposure1_lo_ci <- res_covidB2[1] - 1.96 *res_covidB2[3]
exposure2_lo_ci <- res_covidB2[2] - 1.96 *res_covidB2[4]
exposure1_up_ci <- res_covidB2[1] + 1.96 *res_covidB2[3]
exposure2_up_ci <- res_covidB2[2] + 1.96 *res_covidB2[4]
exposure1_or = exp(res_covidB2[1])
exposure2_or = exp(res_covidB2[2])
exposure1_or_lci95 = exp(exposure1_lo_ci)
exposure2_or_lci95 = exp(exposure2_lo_ci)
exposure1_or_uci95 = exp(exposure1_up_ci)
exposure2_or_uci95 = exp(exposure2_up_ci)

res_covidB2_OR <- cbind(exposure1_or,exposure2_or,exposure1_or_lci95,exposure2_or_lci95,exposure1_or_uci95,exposure2_or_uci95)
#Multivariable MR

#Estimate Std. Error  t value     Pr(>|t|)
#exposure1 0.31996108 0.04712511 6.789610 2.712429e-11
#exposure2 0.02896687 0.02237125 1.294826 1.958793e-01


# ========================================================================================== ## 
##                              Robust causal effect estimation

robust_res_covidB2<- qhet_mvmr(F.data, mvmrcovmatrix, CI = F, iterations = 100)


# ====================================================================================================================== ## 
# ====================================================================================================================== ## 
# ====================================================================================================================== ## 
# ====================================================================================================================== ##

##                                      COVID-B2
# ========================================================================================== ## 
##                                      OUTCOME
##  Read in covidA2 outcome Files
covidC2_ss <- read_tsv(covidC2_ss.path)%>% 
  select(SNP,BETA,SE,TRAIT, P ) 

covidc2_bmi_t2d_snps <- covidC2_ss %>% filter(SNP %in% snps)%>% 
  rename(COVIDC2_beta = BETA, COVIDC2_se = SE)

covidc2_mvmr <-full_join(exposure_mvmr,covidc2_bmi_t2d_snps, by = "SNP" )%>%
  select(BMI_beta,BMI_se,T2D_beta, T2D_se, COVIDC2_beta ,COVIDC2_se, SNP ) 

covidc2_mvmr<- na.omit(covidc2_mvmr)


# ========================================================================================== ## 
##                                  Format summary data     

covidc2_mvmr.data <- format_mvmr(BXGs = covidc2_mvmr[,c(1,3)],
                      BYG = covidc2_mvmr[,5],
                      seBXGs = covidc2_mvmr[,c(2,4)],
                      seBYG = covidc2_mvmr[,6],
                      RSID = covidc2_mvmr[,7])

# ========================================================================================== ## 
##                                  Test for weak instruments 

sres <- strength_mvmr(r_input = covidc2_mvmr.data, gencov = 0)


#Conditional F-statistics for instrument strength

#exposure1 exposure2
#F-statistic  30.57493  15.61036

mvmrcovmatrix<-matrix(c(1,-0.1,-0.05,-0.1,1,0.2,-0.05,0.2,1), nrow = 3, ncol = 3)
Xcovmat<-phenocov_mvmr(mvmrcovmatrix,covidc2_mvmr.data[,5:7])
sres2 <- strength_mvmr(r_input = covidc2_mvmr.data, gencov = Xcovmat)


# ========================================================================================== ## 
##            Test for horizontal pleiotropy using conventional Q-statistic estimation

pres <- pleiotropy_mvmr(r_input = covidc2_mvmr.data, gencov = 0)
# Q-Statistic for instrument validity:


# ========================================================================================== ## 
##                              Estimate causal effects

res_covidC2 <- ivw_mvmr(r_input = covidc2_mvmr.data)
exposure1_lo_ci <- res_covidC2[1] - 1.96 *res_covidC2[3]
exposure2_lo_ci <- res_covidC2[2] - 1.96 *res_covidC2[4]
exposure1_up_ci <- res_covidC2[1] + 1.96 *res_covidC2[3]
exposure2_up_ci <- res_covidC2[2] + 1.96 *res_covidC2[4]
exposure1_or = exp(res_covidC2[1])
exposure2_or = exp(res_covidC2[2])
exposure1_or_lci95 = exp(exposure1_lo_ci)
exposure2_or_lci95 = exp(exposure2_lo_ci)
exposure1_or_uci95 = exp(exposure1_up_ci)
exposure2_or_uci95 = exp(exposure2_up_ci)

res_covidC2_OR <- cbind(exposure1_or,exposure2_or,exposure1_or_lci95,exposure2_or_lci95,exposure1_or_uci95,exposure2_or_uci95)


#Multivariable MR

#Estimate  Std. Error  t value     Pr(>|t|)
#exposure1 0.106996292 0.018683544 5.726766 1.619461e-08
#exposure2 0.008812638 0.008916722 0.988327 3.233909e-01
# ========================================================================================== ## 
##                              Robust causal effect estimation

robust_res_covidC2<- qhet_mvmr(F.data, mvmrcovmatrix, CI = F, iterations = 100)






