library(here)
library(tidyverse)

here::i_am("code/dataprep_pooled.R")
output.path <- here("output/dataprep_pooled")
if (!dir.exists(here(output.path))) dir.create(here(output.path))


pooled.meta <- read_csv(here('output/table2/pooled_cohort.csv'))
# MAF DATA PREP ####


# JHU MAF
jhu.maf <- read_tsv(here('output/dataprep_jhu/jhu_cohort.maf'))
jhu.meta <- read_csv(here('output/supfig1/suptable1.csv'))
jhu.meta <- jhu.meta %>% filter(!is.na(ancestry_group))
jhu.maf <- left_join(jhu.maf, dplyr::select(jhu.meta, ptID, sampleID), by = c("Tumor_Sample_Barcode" = "sampleID"))

# TCGA MAF 
tcga.maf <- read_tsv(here('output/dataprep_TCGA/matched_tcga_maf.tsv'))
tcga.maf <- tcga.maf %>% mutate(ptID = sub("^(([^-]*-){2}[^-]*).*", "\\1", Tumor_Sample_Barcode),
                                sampleID = sub("^(([^-]*-){3}[^-]*).*", "\\1", Tumor_Sample_Barcode))


## Bind Mutation Data
overlap.maf <- intersect(colnames(jhu.maf), colnames(tcga.maf))
maf <- rbind(select(jhu.maf, all_of(overlap.maf)), select(tcga.maf, all_of(overlap.maf)))
maf %>% write_tsv(here(output.path, 'pooled_maf.tsv'))

# CNA DATA PREP ####
jhu.cns <- read_rds(here('output/dataprep_jhu/cna_cytoband.rds'))
tcga.cns <- readRDS(here('output/dataprep_TCGA/cna_tcga_cytoband.rds')) 

hist(assays(tcga.cns["chr3p21.1", ])[[1]], xlim = c(-1,0.5))
hist(assays(jhu.cns["chr3p21.1", ])[[1]], xlim = c(-1,0.5))

hist(assays(tcga.cns["chr5q35.3", ])[[1]], xlim = c(-1,1))
hist(assays(jhu.cns["chr5q35.3", ])[[1]], xlim = c(-1,1))

# Pipeline recommended cuttoff of -0.3 for del and 0.3 for gain is reasonable in both datasets 
jhu.cns.mat <- assays(jhu.cns)[[1]]
tcga.cns.mat <- assays(tcga.cns)[[1]]
colnames(jhu.cns.mat) <- jhu.cns$ptID
colnames(tcga.cns.mat) <- sub("^(([^-]*-){2}[^-]*).*", "\\1", colnames(tcga.cns.mat))

cns.mat <- cbind(jhu.cns.mat, tcga.cns.mat)
pooled.cns.meta <- pooled.meta %>% filter(ptID %in% colnames(cns.mat))

cns.mat <- cns.mat[,pooled.cns.meta$ptID]
pooled.cns.meta <- pooled.cns.meta %>% dplyr::select(ptID, everything())

cns.se <- SummarizedExperiment(assays = list(cns = cns.mat), colData = pooled.cns.meta)

write_rds(cns.se, here(output.path, "pooled_cns_se.rds"))

# RNAseq DATA PREP ####
# Batch correction to be performed


## Loading JHU Cohort
jhu.se <- readRDS(here("output/dataprep_jhu/jhu_se.rds"))
keep <- rowSums(assays(jhu.se)[[1]]) != 0
jhu.se <- jhu.se[keep, jhu.se$Type == "Tumor"]

## Loading TCGA Cohort
tcga.se <- readRDS(here("output/dataprep_TCGA/gene_expression_tcga_se.rds"))
keep <- rowSums(assays(tcga.se)[[1]]) != 0
tcga.se <- tcga.se[keep,tcga.se$tissue_type == "Tumor"]

## Bind on overlapping genes
overlap <- intersect(rownames(assays(jhu.se)[[1]]), rownames(assays(tcga.se)[[1]]))


jhu.se <- jhu.se[overlap,] 
tcga.se <- tcga.se[overlap,] 

jhu.counts <- assays(jhu.se)[["counts"]]
colnames(jhu.counts) <- jhu.se$ptID
jhu.tpm <- assays(jhu.se)[["tpm"]]
colnames(jhu.tpm) <- jhu.se$ptID

tcga.counts <- assays(tcga.se)[["counts"]]
colnames(tcga.counts) <- sub("^(([^-]*-){2}[^-]*).*", "\\1", colnames(tcga.counts))

tcga.tpm <- assays(tcga.se)[["tpm"]]
colnames(tcga.tpm) <- sub("^(([^-]*-){2}[^-]*).*", "\\1", colnames(tcga.tpm))

pooled.tpm <- cbind(jhu.tpm, tcga.tpm)
pooled.counts <- cbind(jhu.counts, tcga.counts)

rna.meta <- pooled.meta %>% filter(ptID %in% colnames(pooled.tpm))

## Batch Correction ####
library(sva)

log2tpm <- log2(pooled.tpm + 1)
rna.meta$cohort <- factor(rna.meta$cohort)
bc.log2.tpm <- ComBat(dat = log2tpm, batch = rna.meta$cohort, par.prior = TRUE)

se <- SummarizedExperiment(assays = list(counts = pooled.counts, tpm = pooled.tpm, "bc.log2.tpm" = bc.log2.tpm), colData = rna.meta)
write_rds(se, here(output.path, "bc_log2tpm_pooled_se.rds"))










