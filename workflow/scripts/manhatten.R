library(tidyverse)
library(TwoSampleMR)
library(ggman)

## Combined manhatten plot

bv.raw <- read_tsv('data/GWAS/Elliott2018bv009.txt.gz')
icv.raw <- read_tsv('data/GWAS/Adams2014icv_fixed.gz')
hc.raw <- read_table2('data/GWAS/Haworth2019hc.chrall.CPRA_b37.tsv.gz', comment = '##')
hcicv.raw <- read_table2('data/GWAS/Haworth2019hcicv.chrall.CPRA_b37.tsv.gz', comment = '##')
hcicvbc.raw <- read_table2('data/GWAS/andrews2019hcicvbc.txt.gz')

bv.dat <- hcicvbc.raw %>% 
  filter(str_detect(SNP, 'rs')) %>%
  select(SNP, CHROM, POS, P) %>% 
  full_join(select(bv.raw, RSID, CHROM, POS, P), by = c('SNP' = 'RSID', 'CHROM', 'POS'), suffix = c('.ICV_HC_BC', '.BV')) %>% 
  full_join(select(hc.raw, DBSNP_ID, CHROM, POS, P), by = c('SNP' = 'DBSNP_ID', 'CHROM', 'POS')) %>% 
  full_join(select(hcicv.raw, DBSNP_ID, CHROM, POS, P), by = c('SNP' = 'DBSNP_ID', 'CHROM', 'POS'), suffix = c('.HC', '.HC_ICV')) %>%
  full_join(filter(icv.raw, str_detect(SNP, 'rs')) %>% select(SNP, CHR, POS, P), by = c('SNP', 'CHROM' = 'CHR', 'POS')) %>% 
  rename(P.ICV = P)


## ============================================= ##
##        Calculate cumulitive position
## ============================================= ##

message("Calculating Positions...")

don <- bv.dat %>%
  
  # Compute chromosome size
  group_by(CHROM) %>%
  summarise(chr_len = max(POS)) %>%
  
  # Calculate cumulative POS of each CHR
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(bv.dat, ., by = c("CHROM" = "CHROM")) %>%
  
  # Add a cumulative POS of each SNP
  arrange(CHROM, POS) %>%
  #mutate(BPcum = 1:nrow(dat)) ## No gaps in chromosome
  mutate(BPcum = POS + tot) ## Gaps in chromsome

axisdf = don %>%
  group_by(CHROM) %>%
  dplyr::summarize(center = (max(BPcum) + min(BPcum)) / 2)

dat.p <- don %>%
  pivot_longer(c('P.ICV_HC_BC', 'P.BV', 'P.HC', 'P.HC_ICV', 'P.ICV'), names_to = 'Trait', values_to = 'P') %>%
  arrange(BPcum, Trait) %>%
  mutate(chr.col = ifelse(CHROM %% 2 == 1, "odd", "even"), 
         Trait = str_replace_all(Trait, 'P.', ""), 
         Trait = fct_relevel(Trait, 'BV', 'HC', 'ICV', 'HC_ICV', 'ICV_HC_BC')) %>%
  filter(!is.na(P))

## ============================================= ##
##        Lead SNPs
## ============================================= ##

bv.lead <- bv.raw %>% 
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
              phenotype_col = 'TRAIT', 
              chr_col = 'CHROM', 
              pos_col = 'POS') %>%
  clump_data(.) %>% 
  as_tibble() %>% 
  left_join(select(don, SNP, tot, BPcum)) %>% 
  select(SNP, CHR = chr.exposure, POS =pos.exposure, P = pval.exposure, tot, BPcum, Trait = exposure)

icv.lead <- icv.raw %>% 
  filter(P < 5e-8) %>% 
  mutate(TRAIT = 'ICV') %>%
  format_data(., type = 'exposure',
              snp_col = 'SNP',
              beta_col = "Beta",
              se_col = "SE",
              eaf_col = "EAF",
              effect_allele_col = "EA",
              other_allele_col = "OA",
              pval_col = "P",
              z_col = "Z",
              samplesize_col = "N",
              phenotype_col = 'TRAIT', 
              chr_col = 'CHR', 
              pos_col = 'POS') %>%
  clump_data(.) %>% 
  as_tibble() %>% 
  left_join(select(don, SNP, tot, BPcum)) %>% 
  select(SNP, CHR = chr.exposure, POS = pos.exposure, P = pval.exposure, tot, BPcum, Trait = exposure)

hc.lead <- hc.raw %>% 
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
              pos_col = 'POS') %>%
  clump_data(.) %>% 
  as_tibble() %>% 
  left_join(select(don, SNP, tot, BPcum)) %>% 
  select(SNP, CHR = chr.exposure, POS =pos.exposure, P = pval.exposure, tot, BPcum, Trait = exposure)

hcicv.lead <- hcicv.raw %>% 
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
              pos_col = 'POS') %>%
  clump_data(.) %>% 
  as_tibble() %>% 
  left_join(select(don, SNP, tot, BPcum)) %>% 
  select(SNP, CHR = chr.exposure, POS =pos.exposure, P = pval.exposure, tot, BPcum, Trait = exposure)

hcicvbc.lead <- hcicvbc.raw %>% 
  filter(P < 5e-8) %>% 
  mutate(TRAIT = 'ICV_HC_BC') %>% 
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
  as_tibble() %>% 
  left_join(select(don, SNP, tot, BPcum)) %>% 
  select(SNP, CHR = chr.exposure, POS =pos.exposure, P = pval.exposure, tot, BPcum, Trait = exposure)

lead.snps <- bind_rows(bv.lead, icv.lead, hc.lead, hcicv.lead, hcicvbc.lead) %>% 
  mutate(Trait = fct_relevel(Trait, 'BV', 'HC', 'ICV', 'HC_ICV', 'ICV_HC_BC'))

## Sample 1% of SNPs with P > 0.05
pval = 0.05
sample_snps <- . %>% 
  filter(P > pval) %>%
  group_by(CHROM, Trait) %>% 
  sample_frac(0.01) %>% 
  ungroup()

p2 <- ggplot() +
  # Horiztonal lines for pvalues
  geom_hline(yintercept = -log10(5e-8), colour = "red", linetype = 3, size = 0.5) +
  geom_hline(yintercept = c(5, 10, 15, 20, 25), colour = 'grey90', size = 0.5) + 
  # background manhatten - p < pval
  geom_point(data = filter(dat.p, P < pval),
             aes(x=BPcum, y=-log10(P), colour = chr.col), size = 0.25) +
  # # background manhatten - p > pval
  geom_point(data = dat.p %>% sample_snps,
             aes(x=BPcum, y=-log10(P), colour = chr.col), size = 0.25) +
  # Highlights for significant loci
  # geom_point(data=highlights, aes(x=BPcum, y=-log10(P), colour = Study),
  #            size=0.5) +
  # Points for lead SNP in each loci
  geom_point(data = lead.snps, aes(x=BPcum, y=-log10(P), colour = Trait), 
             size=2.5, shape = 18) +
  facet_wrap(. ~ Trait, ncol = 1) + 
  # # Locus Labels
  # geom_text(data = new_lab, aes(x=x, y=y_lab, label = loci_new),
  #           angle = 90, hjust = 1, size = 3) +
  # # Lines between locus label and locus
  # geom_segment(data = new_lab, 
  #              aes(x = BPcum, xend = x, y = -1, yend = y_seg), size = 0.1) +
  # # Points for if locus significant in Study
  # geom_point(data = dat_p, 
  #            aes(x = x, y = study.y, colour = Study, shape = present), size=2.5) +
  # # Labels for Study 
  # geom_text(data = study_lab , aes(x=x, y=study.y, label = Study2),
  #           hjust = 1, size = 3, nudge_x = -30000000) + 
  geom_hline(yintercept = 0, colour = 'black', size = 0.5) +
  scale_x_continuous(label = axisdf$CHROM,
                     breaks= axisdf$center ) +
  scale_y_continuous(breaks = c(0, 10, 20)) +
  scale_color_manual(values = c("odd" = "grey50", "even" = "grey25",
                                "BV" = "#377EB8",
                                "HC" = "#4DAF4A",
                                "ICV" = "#984EA3",
                                "HC_ICV" = "#FF7F00",
                                "ICV_HC_BC" = "#E41A1C")) +
  theme_bw() +
   labs(y=bquote("-"~log[10]~"("~italic(p)~")")) +
  theme(
    panel.background = element_rect(fill = "white", colour = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position="none",
    panel.border = element_blank(),
    #panel.grid.major.x = element_blank(),
    #panel.grid.minor.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x=element_blank(),
    axis.text.x=element_blank(),
    text = element_text(size=10),
  )  + scale_shape_manual(values=c(19,1))
p2 
ggsave("results/plots/bv_manhattan.jpg", plot = p2, width = 7.5, height = 10, units = "in")


  