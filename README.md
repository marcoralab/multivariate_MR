# MVMR: BMI + T2D -> SARS-CoV-2 phenotypes
Multivariable MR to investigate the direct effects of body mass index (BMI) and type II diabetes (T2D) on the risk of SARS-CoV-2 phenotypes

```/sc/arion/projects/LOAD/harern01/projects/mvmr_bmi_T2D_COVID```

## Data Sources

**Body Mass Index**
- Yengo L., et al. (2018).Meta-analysis of genome-wide association studies for height and body mass index in ∼700000 individuals of European ancestry.[Human Molecular Genetics 27(20),3641–3649] (https://doi.org/10.1093/hmg/ddy271)

**Type II diabetes**
- Mahajan A., et al. (2018). Fine-mapping type 2 diabetes loci to single-variant resolution using high-density imputation and islet-specific epigenome maps.[Nat Genet 50, 1505–1513]. (https://doi.org/10.1038/s41588-018-0241-6)

**SARS-CoV-2 phenotypes (COVID-19 HGI data freeze 6)**
- COVID19_HGI_A2: Very severe respiratory confirmed covid (8,779 samples) vs. population (1,001,875 samples)
- COVID19_HGI_B2: Hospitalized covid (14,480 samples) vs. not hospitalized covid (73,191 samples)
- COVID19_HGI_C2: Reported infection of Covid (112,612 samples) vs. population (2,474,079 samples)

**Study Characteristics**

| Author | trait | n samples | n SNPs | n indep. SNPs |
| ------ | ----- | --------- | ------ | ------------- |
| Yengo | BMI | 690,495 | 2,336,258 | 516 |
| Mahajan | Diabetes | 898,130 | 23,321,292 | 201 |
| COVID-HGI | COVID: A2 | 1,639,838 | 10,071,637 | - |
| COVID-HGI | COVID: B2 | 2,509,514 | 10,035,826 | - |
| COVID-HGI | COVID: C2 | 2,393,659 | 10,563,782 | - |


**Multivariable Mendelian Randomization**

The genetic correlation, heritability estimates, and Mendelian Randomization were performed as per methods described in Niemi et. al 2021. Multivariable Mendelian Randomization (MVMR) was used to estimate the direct effects of body mass index (BMI) and type II diabetes (T2D) on the risk of SARS-CoV-2 phenotypes, by including both exposures within the same model. We selected all independent (r2 = 0.001; kb = 10000) genome-wide significant (p < 5x10-8) SNPs associated with BMI and type II diabetes; and using the full list of SNPs associated with both BMI and T2D, performed a second clumping procedure (r2 = 0.001; kb = 10000) to obtain independent SNPs. After clumping the full list of SNPs from both BMI and T2D and restricting to SNPs found in the SARS-CoV-2 GWASs, a total of 721 SNPs were available for multivariable MR (516 were associated with BMI only, 201 were associated with T2D only, and 4 SNPs overlap between both GWAS). MVMR was performed using the “MVMR” package (Sanderson et. al 2019). The sensitivity analyses conducted using “MVMR” require estimates of the pairwise covariances between each instrument and each exposure, as such, we used the phenotypic correlation between BMI and T2D and summary data to generate estimates of the covariances. Phenotypic correlations between BMI and T2D were estimated from LDSC regression intercept (r2 = 0.13)(Zheng et. al 2018). Next, we calculated conditional F-statistics to evaluate the presence of weak instruments. The conditional F-statistic for both BMI and T2D were >10, indicating that the selected instruments were strongly associated with their corresponding exposure (Sanderson et. al 2019,Sanderson et. al 2021). We then assessed multivariable instrument pleiotropy using the modified Cochran’s Q-statistic that accounts for potential weak instrument bias, with evidence of heterogeneity indicative of a violation in the escalation restriction assumption in MR (Sanderson et. al 2019,Sanderson et. al 2021). For each of the SARS-CoV-2 phenotypes MVMR models there was evidence of heterogeneity, suggesting the causal estimates from IVW-MVMR may be biased due to the presence of horizontal pleiotropy (Table S10c). Inverse weighted multivariable MR model was used to estimate the direct effect of BMI and T2D upon each of the SARS-CoV-2 phenotypes. When MVMR assumptions are violated, as indicated by the presence of weak instruments or horizontal pleiotropy, it is possible to obtain more robust causal estimates Q- statistic to minimization.
