library(here)
library(SummarizedExperiment)
library(tidyverse)
here::i_am("code/dataprep_TCGA.R")
output.path <- here("output/dataprep_TCGA")
if (!dir.exists(here(output.path))) dir.create(here(output.path))

source(here('code/functions/fileprep_cns_to_binned_matrix.R'))

# PREPARE CNS DATA ------------------------------------------------------------

## Load Propensity Matched Data ####
meta <- read_tsv(here('output/supfig2/matched_kirc_meta.tsv'))

## Load CNS data ####
cns <- read_rds(here('data/tcga/KIRC_DNAcopy_segment_CN.rds'))
cns <- cns %>% dplyr::rename(caseID = Sample)

### Clean/Reformt CNS Manifest
cns.meta <- read_csv('data/tcga/KIRC_DNAcopy_segment_CN_manifest.csv')
cns.meta <- cns.meta %>% dplyr::rename(ptID = cases.submitter_id, sampleID = sample.submitter_id, caseID = cases)

### Add sampleID to CNS Manifest 
cns <- left_join(cns, select(cns.meta, sampleID, caseID))

### Add CNS metadata to  Propensity matched cohort
meta <- left_join(meta, cns.meta)

### Remove duplicate files 
meta <- meta %>% group_by(sampleID) %>% slice_head()


### Exclude additional Primary 
meta <- meta %>% filter(sample_type != "Additional - New Primary")
meta <- meta %>% mutate(Type = if_else(str_detect(sample_type, "Normal"), "Normal", "Tumor"))
meta %>% ungroup %>% dplyr::count(Type)

### Exclude normal samples
meta <- meta %>% filter(Type != "Normal")

### Subset CNS data for matched cohort
cns <- cns %>% filter(sampleID %in% meta$sampleID)


## Create SE Object ####

### Reformat Chromosome variable for compatibility with cytobands
cns <- cns %>% mutate(Chromosome = paste0("chr",Chromosome))

cytoband <- read_tsv(here("data/dna/cytoband_table.tsv"))

### Generate binned cytoband object
cytoband <- cytoband %>% dplyr::rename(chr = `chrom`, start = chromStart, end = chromEnd)
cytoband$cytoband <- str_c(cytoband$chr, cytoband$name)

### Subset problematic centromeric regions
cytoband <- cytoband %>% filter(chr %in% paste0("chr", 1:22), !gieStain %in% c("acen", "gvar"))
cns <- cns %>% filter(Chromosome %in%  paste0("chr", 1:22))

### granges object and format to mat
cytoband_gr = GRanges(seqnames = cytoband$chr, ranges = IRanges(cytoband$start + 1, cytoband$end), mcols = cytoband)
mat <- fileprep_cns_to_binned_matrix(cns.object = cns, gr_object = cytoband_gr, grange.mcol = "mcols.cytoband", value = "Segment_Mean")

### perpare summarized experiment object
meta$sampleID <- factor(meta$sampleID, levels = colnames(mat))
meta <- meta %>% arrange(sampleID )
meta <- meta %>% select(sampleID, everything())

cytoband <- SummarizedExperiment(assays = list("cytoband_log2fold" = mat), colData = as.data.frame(meta), rowData = cytoband) 
rowData(cytoband)$cytoband <- factor(rowData(cytoband)$cytoband)

###  Save 
cytoband %>% saveRDS(here(output.path,"cna_tcga_cytoband.rds"))



# Prepare Gene Expression Data ---------------------------------------------------

## Load Propensity Matched Data
meta <- read_tsv(here('output/supfig2/matched_kirc_meta.tsv'))

## Load RNA data
se.rna <- read_rds(here('data/tcga/kirc_rnaseq_se.rds'))
rna.meta <- as.data.frame(colData(se.rna))
rna.meta <- rna.meta %>% dplyr::rename(ptID = patient, caseID = barcode, sampleID = sample)

##Subset for GeneNames
se.rowdata <- as.data.frame(rowData(se.rna))
se.rowdata <- se.rowdata %>% filter(!is.na(gene_name)) %>% dplyr::group_by(gene_name) %>% slice_head()



## Pull counts and TPM
counts <- assays(se.rna)[["unstranded"]]
tpm <- assays(se.rna)[["tpm_unstrand"]]


## Subset for selected gene names
se.rowdata <- as.data.frame(rowData(se.rna))
se.rowdata <- se.rowdata %>% filter(!is.na(gene_name)) %>% dplyr::group_by(gene_name) %>% slice_head()
counts <- counts[se.rowdata$gene_id,]
tpm <- tpm[se.rowdata$gene_id,]

##Rename rows
rownames(counts) <- se.rowdata$gene_name
rownames(tpm) <- se.rowdata$gene_name

## Subset for Propensity = matched cohort
rna.meta <- rna.meta %>% filter(ptID %in% meta$ptID, shortLetterCode != "TAP")
rna.meta %>% dplyr::count(shortLetterCode) 
rna.meta <- rna.meta %>% mutate(Type = if_else(shortLetterCode == "TP", "Tumor", "Normal"))

rna.meta <- rna.meta %>% group_by(sampleID) %>% slice_head() 
rna.meta <- left_join(rna.meta, meta)

counts <- counts[,rna.meta$caseID]
tpm <- tpm[,rna.meta$caseID]

colnames(counts) <-  rna.meta$sampleID
colnames(tpm) <-  rna.meta$sampleID



## Create SE object ####
se <- SummarizedExperiment(assays = list("counts" = counts, "tpm" = tpm), colData = as.data.frame(rna.meta), rowData = se.rowdata) 
se %>% saveRDS(here(output.path,"gene_expression_tcga_se.rds"))


# Prepare SNV Data ####
maf <- read_rds(here('data/tcga/kirc_maf.rds'))
meta <- read_tsv(here('output/supfig2/matched_kirc_meta.tsv'))

# Add sample and ptID
maf <- maf %>%
  mutate(
    ptID = sub("^(([^-]*-){2}[^-]*).*", "\\1", Tumor_Sample_Barcode),
    sampleID = sub("^(([^-]*-){3}[^-]*).*", "\\1", Tumor_Sample_Barcode),
    analyteID = str_sub(Tumor_Sample_Barcode, 20, 20)
  )

maf %>% dplyr::count(analyteID)
maf <- maf %>% filter(analyteID == "D", ptID %in% meta$ptID)
meta <- left_join(meta, unique(select(maf, ptID, Tumor_Sample_Barcode)))
meta <- meta %>% group_by(ptID) %>% slice_head()
maf <- maf %>% filter(Tumor_Sample_Barcode %in% meta$Tumor_Sample_Barcode)

#' Subset for DNA analyte rather than whole genome amplification [amplification](https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/portion-analyte-codes)
maf %>% write_tsv(here(output.path, "matched_tcga_maf.tsv"))

sessionInfo()








