library(here)
library(tidyverse)
library(readxl)


here::i_am("code/dataprep_jhu.R")
output.path <- here("output/dataprep_jhu")
if (!dir.exists(here(output.path))) dir.create(here(output.path))


#Prepare JHU MAF, Copy Number, and RNAseq Data ####

## MAF ####
meta <- read_csv(here('output/supfig1/suptable1.csv'))
id.key <- read_excel(here('data/sample_annotation/ID_Key.xlsx'))

data.path <- here("data/dna/mutect2_maf/")
files <- list.files(here(data.path))


output <- list()
ncols <- vector()

for (i in seq_along(1:length(files))){
  maf <- read_tsv(here(data.path, files[i]), col_types = cols(.default = "c"))
  output[[i]] <- maf
  ncols[i] <- ncol(maf)
}

maf <- do.call(bind_rows, output) 

maf <- left_join(maf, select(id.key, Tumor_Sample_Barcode, sampleID))
maf <- maf %>% select(-Tumor_Sample_Barcode)
maf <- maf %>% dplyr::rename(Tumor_Sample_Barcode = sampleID)

maf <- maf %>% mutate(Matched_Norm_Sample_Barcode = if_else(Matched_Norm_Sample_Barcode == "NORMAL", "Pooled Normal", "Matched Normal"))

maf <- maf %>% filter(Tumor_Sample_Barcode %in% meta$sampleID)


# Filter Pooled Normal Pipeline ####
pooled.normal.maf <- maf %>% filter(Matched_Norm_Sample_Barcode == "Pooled Normal")


#' First Filter:  Allelic Frequency < 0.01 (1%) in 1000 Genomes AND gnomAD OR not detected in reference databases
#' Results in 146600 variants
 
pooled.normal.maf <- pooled.normal.maf %>%
  filter((AF < 0.01 | is.na(AF)) & (gnomAD_AF < 0.01 | is.na(gnomAD_AF))) 
#' Next, Filter for coding domain sequences defined below: 
#' Results in 8527 variants
cds <-
  c("Frame_Shift_Del",
    "Frame_Shift_Ins",
    "In_Frame_Del",
    "In_Frame_Ins",
    "Missense_Mutation",
    "Nonsense_Mutation",
    "Nonstop_Mutation",
    "Translation_Start_Site"
  )

pooled.normal.maf <- pooled.normal.maf %>%
  filter(Variant_Classification %in% cds)

#' Next, subset missense mutations and filter only for likely deleterious mutations
#' #This process takes 5104 missense mutations and reduces to 2606
#' The total number of Variants after recombining is 6029
missense <- pooled.normal.maf %>% filter(Variant_Classification == "Missense_Mutation")
missense <- missense %>% filter(str_detect(SIFT, "deleterious") | str_detect(CLIN_SIG,"pathogenic|uncertain|conflicting") | str_detect(PolyPhen, "damaging")) 

pooled.normal.maf <- pooled.normal.maf %>% filter(Variant_Classification != "Missense_Mutation")
pooled.normal.maf <- bind_rows(pooled.normal.maf, missense)

pooled.normal.maf %>% dplyr::count(Tumor_Sample_Barcode)

#' Next Filter the Matched.maf
#' Starts with 22006 Variants
#' Filtering for coding domain mutations results in 8497 Variants
matched.maf <- maf %>% filter(Matched_Norm_Sample_Barcode != "Pooled Normal")

cds <-
  c("Frame_Shift_Del",
    "Frame_Shift_Ins",
    "In_Frame_Del",
    "In_Frame_Ins",
    "Missense_Mutation",
    "Nonsense_Mutation",
    "Nonstop_Mutation",
    "Translation_Start_Site",
    "Splice_Site"
  )

matched.maf <- matched.maf %>% filter(Variant_Classification %in% cds)

maf.df <- bind_rows(matched.maf, pooled.normal.maf)

# Subset for samples which passed QC ####
wes_qc <- read_excel(here('data/dna/wes_qc/QC.xlsx'))
wes_qc <- wes_qc %>% filter(WES_STATUS == "PASSED")
maf.df <- maf.df %>% filter(Tumor_Sample_Barcode %in% wes_qc$sampleID)
write_tsv(maf.df, here(output.path, "jhu_cohort.maf"))


# Copy Number ####


data.path <- here("data/dna/cnvkit_cns/")
files <- list.files(path = data.path, pattern = ".cns", recursive = "TRUE")

data_list <- list()
for(i in seq_along(files)){
  sample_cns <- read_tsv(paste0(data.path,files[i]))
  vector <- strsplit(files[i], "_",fixed=TRUE)[[1]][1:2]
  sample_id <- str_c(vector, collapse = "_")
  sample_cns$Tumor_Sample_Barcode <- sample_id
  data_list[[i]] <- sample_cns
  data_list
}

cns <- do.call(rbind, data_list)
cns$length <- cns$end-cns$start


# Filter samples which failed quality control
cns <- left_join(cns, select(id.key, Tumor_Sample_Barcode, sampleID))
cns <- cns %>% select(-Tumor_Sample_Barcode)
meta <- meta %>% filter(WES_STATUS == "PASSED", Type == "Tumor")
cns <- cns %>% filter(sampleID %in% meta$sampleID)
cns <- cns %>% dplyr::rename(Chromosome = chromosome, Start = start, End = end)
write_tsv(cns, here(output.path, "jhu_merged.cns"))


source(here('code/functions/fileprep_cns_to_binned_matrix.R'))

# #Bin CNS file
# session <- browserSession()
# query <- rtracklayer::ucscTableQuery("hg38", table = "cytoBandIdeo")
# cytoband_table <- rtracklayer::getTable(query)
# write_tsv(cytoband_table, here("data/dna/cytoband_table.tsv"))

cytoband <- read_tsv(here("data/dna/cytoband_table.tsv"))
cytoband <- cytoband %>% dplyr::rename(chr = `chrom`, start = chromStart, end = chromEnd)
cytoband$cytoband <- str_c(cytoband$chr, cytoband$name)

## Filter for autosomes, and remove noisy centromeric and giestain variable regions

cytoband <- cytoband %>% filter(chr %in% paste0("chr", 1:22), !gieStain %in% c("acen", "gvar"))
cns <- cns %>% filter(Chromosome %in%  paste0("chr", 1:22))

# Make Grange object: 
library(GenomicRanges)
cytoband_gr = GRanges(seqnames = cytoband$chr, ranges = IRanges(cytoband$start + 1, cytoband$end), mcols = cytoband)

mat <- fileprep_cns_to_binned_matrix(cns.object = cns, gr_object = cytoband_gr, grange.mcol = "mcols.cytoband")

## Save RDS
meta$sampleID <- factor(meta$sampleID, levels = colnames(mat))
meta <- meta %>% arrange(sampleID )
meta <- meta %>% select(sampleID, everything())

library(SummarizedExperiment)
cytoband <- SummarizedExperiment(assays = list("cytoband_log2fold" = mat), colData = as.data.frame(meta), rowData = cytoband) 
rowData(cytoband)$cytoband <- factor(rowData(cytoband)$cytoband)

cytoband %>% saveRDS(here(output.path,"cna_cytoband.rds"))


# RNAseq ####
# Load STAR count outputs ####
data.path <- here('data/rna/star_counts/')
file.names <- list.files(here('data/rna/star_counts/'))

count.output <- list()
qc.output <- list()

for(i in 1:length(file.names)){
  file.path <- paste0(data.path,file.names[i])
  sampleID <- strsplit(file.names[i], "_",fixed=TRUE)[[1]][1]
  
  #pull counts
  data <- read.table(file.path,
                     stringsAsFactors = FALSE,
                     sep = "\t",
                     header = FALSE,
                     row.names = NULL,
                     check.names = FALSE,
                     skip = 4) 
  unstranded <- data[,2]
  names(unstranded) <- data[,1]
  count.output[[i]] <- unstranded
  names(count.output)[i] <- sampleID
  
  
  #Pull qc metrics
  read_metrics <- read.table(file.path,
                             stringsAsFactors = FALSE,
                             sep = "\t",
                             header = FALSE,
                             row.names = NULL,
                             check.names = FALSE,
                             nrows = 4)  
  
  df <- data.frame(t(read_metrics$V2))
  colnames(df) <- read_metrics$V1
  df$sampleID <- sampleID
  df$total_mapped <- sum(unstranded)
  qc.output[[i]] <- df
}

qc <- do.call(bind_rows, qc.output)
qc <- qc %>% select(sampleID, everything()) %>% mutate(
  per_unmapped = N_unmapped/(total_mapped + N_unmapped) * 100
)

counts <- do.call(bind_cols, count.output)
rownames(counts) <- names(count.output[[1]]) 
counts <- as.matrix(counts)


#Generate TPM ####
transcript.lengths <- read_tsv(here('data/rna/jhu_gene_annotation/hg38_gene_length.txt'))
glv <- transcript.lengths$Gene_length
names(glv) <- transcript.lengths$Gene_symbol
glv <- glv[rownames(counts)] #Creates a vector of gene lengths in the same order as my matrix 

library(DGEobj.utils)
tpm <- convertCounts(counts, unit = "tpm", geneLength = glv)


keep <- rowSums(tpm) != 0
log2tpm <- log2(tpm[keep,] + 1)
colMedians(log2tpm) == 0


qc$medlog2tpm0 <- colMedians(log2tpm) == 0
qc <- qc %>% mutate(RNAseq_STATUS = 
                      if_else(medlog2tpm0, "FAILED", "PASSED"))



write_csv(qc, here(output.path, "rnaseq_qc.csv"))
passed <- qc %>% filter(RNAseq_STATUS == "PASSED")
counts.passed <- counts[, passed$sampleID]
tpm.passed <- tpm[, passed$sampleID]


##Sum Exp Object: 
meta <- read_csv(here('output/supfig1/sample_lvl.csv'))
rna.se.meta <- meta %>% filter(sampleID %in% colnames(tpm.passed))
rownames(rna.se.meta) <- rna.se.meta$sampleID

tpm.passed <- tpm.passed[,rna.se.meta$sampleID]
counts.passed <- counts.passed[,rna.se.meta$sampleID]
se <- SummarizedExperiment(assays = list("counts" = counts.passed, "tpm" = tpm.passed),
                           colData = rna.se.meta)

saveRDS(se, here(output.path, "jhu_se.rds"))


sessionInfo()


