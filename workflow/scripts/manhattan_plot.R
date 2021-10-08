#!/usr/bin/Rscript

## Load Packages 
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(ggman))

## Read in arguments
args = commandArgs(trailingOnly = TRUE) # Set arguments from the command line
infile_gwas = args[1]
infile_clump = args[2]
outfile_plot = args[3]
PlotTitle = args[4]

message(paste(PlotTitle, '\n'))

## Read in GWAS and Plink Clumped File
trait.gwas <- read_tsv(infile_gwas, comment = '##', guess_max = 15000000)

trait.clump <- suppressMessages(read_table2(infile_clump)) %>% 
  filter(!is.na(CHR)) %>% 
  select(CHR, F, SNP, BP, P, TOTAL, NSIG)

## Index SNPs
IndexSnps1 <- filter(trait.clump, P <= 5e-8) %>% pull(SNP)

## Calculate Maximum Pvalue
max.p <- max(-log10(filter(trait.clump, P > 0)$P)) + 5

## Plot GWAS
message('PLOTTING MANHATTEN PLOT\n')
p.TRAIT <- ggman(filter(trait.gwas, P < 0.05), snp = 'DBSNP_ID', chrom = 'CHROM', bp = 'POS', pvalue = 'P', ymin = 0, ymax = max.p, 
                 title = PlotTitle, sigLine = -log10(5e-8), relative.positions = TRUE) + 
  theme_classic() + theme(text = element_text(size=10)) 

if(length(IndexSnps1) >= 1){
  message('PLOTTING LEAD SNPS\n')
  p.TRAIT <- ggmanHighlight(p.TRAIT, highlight = IndexSnps1, shape = 18, colour = 'red')
}

message('EXPORTING PLOT\n')
ggsave(p.TRAIT, device = 'png', units = 'in', width = 10, height = 5, plot = out.TRAIT)
