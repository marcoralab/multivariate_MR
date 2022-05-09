'''Snakefile for Multivariable Mendelian Randomization'''

import os
from itertools import product
import pandas as pd

shell.prefix('module load plink/1.90 R/3.5.1; ')

configfile: "config.yaml"
ExposureCode = config["ExposureCode"]
OutcomeCode = config["OutcomeCode"]
Pthreshold = config["Pthreshold"]
mvexp1 = config["mvexp1"]
mvexp2 = config["mvexp2"]
REF = '../MRcovid/data/raw/EUR_All_Chr'
r2 = 0.001
kb = 10000
# Filter forbidden wild card combinations
## https://stackoverflow.com/questions/41185567/how-to-use-expand-in-snakemake-when-some-particular-combinations-of-wildcards-ar
def filter_combinator(combinator, blacklist):
    def filtered_combinator(*args, **kwargs):
        for wc_comb in combinator(*args, **kwargs):
            # Use frozenset instead of tuple
            # in order to accomodate
            # unpredictable wildcard order
            if frozenset(wc_comb) not in blacklist:
                yield wc_comb
    return filtered_combinator

forbidden = {frozenset(wc_comb.items()) for wc_comb in config["missing"]}
filtered_product = filter_combinator(product, forbidden)

rule all:
    input:
        expand("input/{exposure}/{outcome}/{exposure}_{pt}_{outcome}_MR_Results.txt",
        filtered_product, exposure = ExposureCode, pt = Pthreshold, outcome = OutcomeCode),
        expand("results/{mvexp1}_{mvexp2}_{OutcomeCode}_sres_IVMMR.txt",mvexp1 = mvexp1, mvexp2 = mvexp2,
        OutcomeCode = OutcomeCode),
        expand("results/{mvexp1}_{mvexp2}_{OutcomeCode}_pres_IVMMR.txt",mvexp1 = mvexp1, mvexp2 = mvexp2,
        OutcomeCode = OutcomeCode),
        expand("results/{mvexp1}_{mvexp2}_{OutcomeCode}_res_IVMMR.txt",mvexp1 = mvexp1, mvexp2 = mvexp2,
        OutcomeCode = OutcomeCode),
        expand("results/{mvexp1}_{mvexp2}_{OutcomeCode}_robust_res_IVMMR.txt",mvexp1 = mvexp1, mvexp2 = mvexp2,
        OutcomeCode = OutcomeCode)

        # expand("input/{mvexp1}/{mvout}/{mvexp1}_{mvexp2}_{pt}_{mvout}_MVMRresults.txt",
        # mvexp1 = mvexp1, mvexp2 = mvexp2, mvout = mvout, pt = Pthreshold)
        # expand('input/{exposure}/{outcome}/{outcome}.txt', outcome = TEST)

## Extract SNPs to be used as instruments in exposure
rule ExposureSnps:
    input:
        summary = "../MRcovid/data/formated/{exposure}/{exposure}_formated.txt.gz",
        ExposureClump = "../MRcovid/data/formated/{exposure}/{exposure}.clumped"
    output: out = "input/{exposure}/{exposure}_{pt}_SNPs.txt"
    params: Pthreshold = '{pt}'
    conda: "workflow/envs/r.yaml"
    script:
        'workflow/scripts/ExposureData.R'


## Extract exposure instruments from outcome gwas
rule OutcomeSnps:
    input:
        ExposureSummary = rules.ExposureSnps.output.out,
        OutcomeSummary = '../MRcovid/data/formated/{outcome}/{outcome}_formated.txt.gz'
    output:
        snps = "input/{exposure}/{outcome}/{exposure}_{pt}_{outcome}.SNPs.txt",
        missing = "input/{exposure}/{outcome}/{exposure}_{pt}_{outcome}.MissingSNPs.txt",
    conda: "workflow/envs/r.yaml"
    script:
        'workflow/scripts/OutcomeData.R'


## Use plink to identify proxy snps instruments that were not avaliable in the outcome
rule FindProxySnps:
    input: rules.OutcomeSnps.output.missing,
    output: "input/{exposure}/{outcome}/{exposure}_{pt}_{outcome}_Proxys.ld",
    params:
        Outcome = "input/{exposure}/{outcome}/{exposure}_{pt}_{outcome}_Proxys",
        ref = REF
    shell:
        """
        if [ $(wc -l < {input}) -eq 0 ]; then
            touch {output}
          else
           plink --bfile {params.ref} \
           --keep-allele-order \
           --r2 dprime in-phase with-freqs \
           --ld-snp-list {input} \
           --ld-window-r2 0.8 --ld-window-kb 500 --ld-window 1000 --out {params.Outcome}
          fi
"""

## Extract proxy SNPs from outcome gwas
rule ExtractProxySnps:
    input:
        OutcomeSummary = '../MRcovid/data/formated/{outcome}/{outcome}_formated.txt.gz',
        OutcomeProxys = rules.FindProxySnps.output,
        OutcomeSNPs = rules.OutcomeSnps.output.snps,
    output:
        proxies = "input/{exposure}/{outcome}/{exposure}_{pt}_{outcome}_ProxySNPs.txt",
        matched = "input/{exposure}/{outcome}/{exposure}_{pt}_{outcome}_MatchedProxys.csv",
    conda: "workflow/envs/r.yaml"
    script:
        'workflow/scripts/ExtractProxySNPs.R'

## Use TwoSampleMR to harmonize exposure and outcome datasets
rule Harmonize:
    input:
        ExposureSummary = "input/{exposure}/{exposure}_{pt}_SNPs.txt",
        OutcomeSummary = "input/{exposure}/{outcome}/{exposure}_{pt}_{outcome}_ProxySNPs.txt",
        ProxySNPs = "input/{exposure}/{outcome}/{exposure}_{pt}_{outcome}_MatchedProxys.csv"
    output:
        Harmonized = "input/{exposure}/{outcome}/{exposure}_{pt}_{outcome}_MRdat.csv",
    params:
        Pthreshold = '{pt}',
        excposurecode = "{exposure}",
        outcomecode = "{outcome}"
    singularity: "docker://mrcieu/twosamplemr"
    script: 'workflow/scripts/DataHarmonization.R'

rule radialMR:
    input:
        mrdat = rules.Harmonize.output.Harmonized
    output:
        mrdatradial = "input/{exposure}/{outcome}/{exposure}_{pt}_{outcome}_MRdatRadial.csv",
        heterogenity = "input/{exposure}/{outcome}/{exposure}_{pt}_{outcome}_modifiedQ.txt"
    params:
        out = "input/{exposure}/{outcome}/{exposure}_{pt}_{outcome}"
    singularity: "docker://mrcieu/twosamplemr"
    script:
        'workflow/scripts/Radial.R'

## Conduct MR analysis
rule MR_analysis:
    input:
        mrdat = rules.radialMR.output.mrdatradial
    output:
        "input/{exposure}/{outcome}/{exposure}_{pt}_{outcome}_MR_heterogenity.txt",
        "input/{exposure}/{outcome}/{exposure}_{pt}_{outcome}_MR_egger_plei.txt",
        "input/{exposure}/{outcome}/{exposure}_{pt}_{outcome}_MR_Results.txt",
        "input/{exposure}/{outcome}/{exposure}_{pt}_{outcome}_MR_nsnps.txt",
    params:
        out = "input/{exposure}/{outcome}/{exposure}_{pt}_{outcome}"
    singularity: "docker://mrcieu/twosamplemr"
    script:
        'workflow/scripts/MR_analysis.R'

# Conduct MR analysis
rule MVMR_SNPlist:
    input:
        exp1_mrdat = "input/{mvexp1}/{mvout}/{mvexp1}_{pt}_{mvout}_MRdatRadial.csv",
        exp2_mrdat = "input/{mvexp2}/{mvout}/{mvexp2}_{pt}_{mvout}_MRdatRadial.csv",
    output:
        out = "input/{mvexp1}/{mvout}/{mvexp1}_{mvexp2}_{pt}_{mvout}_CombinedSNPlist.txt",
    conda: "envs/r.yaml"
    script: 'workflow/scripts/CombineMVMRExposures.R'

rule MVMRclump:
    input:
        gwas = 'input/{mvexp1}/{mvexp1}.txt',
        mvmr_snplists = "input/{mvexp1}/{mvout}/{mvexp1}_{mvexp2}_{pt}_{mvout}_CombinedSNPlist.txt"
    output:
        exp_clumped = 'input/{mvexp1}/{mvout}/{mvexp1}_{mvexp2}_{pt}_{mvout}.clumped',
    params:
        ref = REF,
        out =  'input/{mvexp1}/{mvout}/{mvexp1}_{mvexp2}_{pt}_{mvout}',
        r2 = r2,
        kb = kb
    shell:
        """
        plink --bfile {params.ref} --keep-allele-order --allow-no-sex  \
        --extract {input.mvmr_snplists} \
        --clump {input.gwas} --clump-snp-field SNP --clump-field P \
        --clump-r2 {params.r2} --clump-kb {params.kb} --clump-p1 1 --clump-p2 1 --out {params.out};
        """

rule MVMR:
    input:
        gwas_exp1 = '../MRcovid/data/formated/{mvexp1}/{mvexp1}_formated.txt.gz',
        gwas_exp2 = '../MRcovid/data/formated/{mvexp2}/{mvexp2}_formated.txt.gz',
        gwas_out = '../MRcovid/data/formated/{mvout}/{mvout}_formated.txt.gz',
        mvmr_snplists = "input/{mvexp1}/{mvout}/{mvexp1}_{mvexp2}_{pt}_{mvout}_CombinedSNPlist.txt"
    output:
        MVMRresults = 'input/{mvexp1}/{mvout}/{mvexp1}_{mvexp2}_{pt}_{mvout}_MVMRresults.txt',
    params:
        mvexp1 = '{mvexp1}',
        mvexp2 = '{mvexp2}',
        mvout = '{mvout}'
    singularity: "docker://mrcieu/twosamplemr"
    script:'workflow/scripts/MVMR2.R'

rule IVW_MVMR:
    input:
        bmi = '../MRcovid/data/formated/{mvexp1}/{mvexp1}_formated.txt.gz',
        t2d = '../MRcovid/data/formated/{mvexp2}/{mvexp2}_formated.txt.gz',
        covid = '../MRcovid/data/formated/{OutcomeCode}/{OutcomeCode}_formated.txt.gz'
    output:
        outfile1 = "results/{mvexp1}_{mvexp2}_{OutcomeCode}_sres_IVMMR.txt",
        outfile2 ="results/{mvexp1}_{mvexp2}_{OutcomeCode}_pres_IVMMR.txt",
        outfile3 ="results/{mvexp1}_{mvexp2}_{OutcomeCode}_res_IVMMR.txt",
        outfile4 ="results/{mvexp1}_{mvexp2}_{OutcomeCode}_robust_res_IVMMR.txt"
    params:
        bmi_snplists = 'input/{mvexp1}/{mvexp1}_5e-8_SNPs.txt',
        t2d_snplists = 'input/{mvexp2}/{mvexp2}_5e-8_SNPs.txt'
    conda: "workflow/envs/r.yaml"
    script:'workflow/scripts/MVMR_IVW.R'
