library(here)
library(tidyverse)
library(readxl)
library(matrixStats)

source(here('code/functions/ggplot_samplelvl_boxplot.R'))

here::i_am("code/supfig1.R")
output.path <- here("output/supfig1")
if (!dir.exists(here(output.path))) dir.create(here(output.path))

#  SUP FIG 1A Flow and numbers ####
sample_data <- read_excel(here('data/sample_annotation/sample_data.xlsx'))
wes_qc <- read_excel(here('data/dna/wes_qc/QC.xlsx'))
rna_qc <- read_csv(here('output/dataprep_jhu/rnaseq_qc.csv'))
ancestry <- read_csv(here('data/dna/admixture/ancestry_estimates.csv'))

#Initial Cohort N 123 samples - N = 66 patients
sample_data %>% dplyr::count(ptID) %>% filter(n == 1) #9 patients with Tumor only
sample_data %>% dplyr::count(ptID) %>% filter(n == 2) #57 patients with matched T/N
sample_data %>% select(ptID, Race) %>% unique() %>% dplyr::count(Race)

## WES PIPELINE ####
wes_qc <- wes_qc %>% filter(sampleID %in% sample_data$sampleID)
wes_qc %>% group_by(Group) %>% dplyr::count(WES_STATUS)

wes_paired <- read_excel(here('data/dna/wes_qc/QC.xlsx'), sheet = 2)
wes_tumor_only <- read_excel(here('data/dna/wes_qc/QC.xlsx'), sheet = 3)

wes_paired <- wes_paired %>% filter(`Paired Analysis` == "PASSED")
wes_tumor_only <- wes_tumor_only %>% filter(`Tumor Only` == "PASSED")

wes_qc <- wes_qc %>% mutate(WES_comparison = case_when(
  sampleID %in% wes_paired$sampleID ~ "Matched Normal",
  sampleID %in% wes_tumor_only$sampleID ~ "Pooled Normal",
  .default = NA
))

wes_qc %>% group_by(Group) %>% dplyr::count(WES_comparison)

# Add WES QC Metrics to Sample Data
wes_qc <- left_join(wes_qc, sample_data)
wes_qc %>% filter(WES_STATUS == "PASSED") %>% group_by(Race) %>% dplyr::count(WES_comparison)
wes_cohort <- wes_qc %>% filter(WES_STATUS == "PASSED", Group == "Tumor")
sample_data <- left_join(sample_data, select(wes_qc, sampleID, Avg_Cov, WES_STATUS, WES_comparison))


## RNAseq Pipeline ####
rna_qc <- rna_qc %>% filter(sampleID %in% sample_data$sampleID)
rna_qc <- rna_qc %>% filter(RNAseq_STATUS == "PASSED")
rna_qc <- left_join(rna_qc, select(sample_data, -year_sugery))

rna_cohort <- rna_qc %>%
  group_by(ptID) %>%
  filter(any(Type == "Tumor")) %>%
  ungroup()


rna_cohort %>% group_by(Race) %>% dplyr::count(Type)

# Add To Sample Data
sample_data <- left_join(sample_data, select(rna_cohort, sampleID, per_unmapped, total_mapped, RNAseq_STATUS))
sample_data <- sample_data %>% mutate(RNAseq_STATUS = if_else(is.na(RNAseq_STATUS), "FAILED", RNAseq_STATUS))

# Add Ancestry Annotation - Note, this is a per patient measurement.
ancestry <- left_join(ancestry, wes_qc)
ancestry <- ancestry %>% group_by(ptID) %>% slice_max(Avg_Cov)
ancestry <- ancestry %>% select(ptID, SAS:AFR)
sample_data <- left_join(sample_data, ancestry)

# Subset for Final Cohort ####
final.cohort <- sample_data %>% filter(!WES_STATUS %in% c("FAILED") | !RNAseq_STATUS %in% c("FAILED"))

#Exclude cases with only passed Normal (7 patients where only normal specimen passed)
final.cohort <- final.cohort %>%
  group_by(ptID) %>%
  filter(any(Type == "Tumor")) %>%
  ungroup()

#Sample Level Statistics
final.cohort %>% dplyr::count(ptID) #107 samples from 59 pts (from initial 123 samples, 66pts)


# Save Patient Level Data ####
pt.lvl <- final.cohort %>% filter(Type == "Tumor")

clinical_data <- read_excel(here('data/clinical_annotation/clinical_annotation.xlsx'))
clinical_data <- clinical_data %>% select(ptID, OS, OS.time, Age, Cr, HTN)

## Add CDK-EPI Score
pt.lvl <- left_join(pt.lvl, clinical_data, by = "ptID")

pt.lvl <- pt.lvl %>%
  mutate(
    k = ifelse(Sex == "Female", 0.7, 0.9),
    alpha = ifelse(Sex == "Female", -0.329, -0.411),
    eGFRcr = 142 * 
      (ifelse(Cr / k < 1, Cr / k, 1)^alpha) *
      (ifelse(Cr / k > 1, Cr / k, 1)^-1.200) * 
      (0.9938^Age) *
      (ifelse(Sex == "Female", 1.012, 1))
  ) %>%
  select(-k, -alpha)

## Add Variable Bins

pt.lvl <- pt.lvl %>% mutate(eGFR_cat = 
                              case_when(eGFRcr < 30 ~ "<30",
                                        eGFRcr >= 30 & eGFRcr <= 45 ~ "30-45",
                                        eGFRcr > 45 & eGFRcr < 60 ~ "45-60",
                                        eGFRcr >= 60 & eGFRcr <= 90~ "61-90",
                                        eGFRcr > 90 ~ ">90"))

pt.lvl <- pt.lvl %>% mutate(GFR_bin = 
                              case_when(eGFRcr <= 60 ~"<60",
                                        eGFRcr > 60 ~">60"))
pt.lvl <- pt.lvl %>% mutate(Age_cat = 
                              case_when(
                                Age < 50 ~ "<50",
                                Age >= 50 & Age <= 60 ~ "50-60",
                                Age > 60 & Age <= 70 ~ "61-70",
                                Age > 70 ~ ">70"
                              ))
pt.lvl <- pt.lvl %>% mutate(HTN = if_else(HTN == 0, "No", "Yes"))

pt.lvl$Age_cat <- factor(pt.lvl$Age_cat, levels = c("<50", "50-60", "61-70", ">70"))
pt.lvl$eGFR_cat <- factor(pt.lvl$eGFR_cat, levels = c("<30", "30-45", "45-60", "61-90", ">90"))

#write_csv(pt.lvl, here(output.path, "pt_lvl_output.csv"))

final.cohort <- left_join(final.cohort, select(pt.lvl, ptID, Age, Cr, HTN, eGFRcr, Age_cat,eGFR_cat, GFR_bin), by = "ptID")

# SupFig1B (WES Average Coverage) ####
wes_qc <- read_excel(here('data/dna/wes_qc/QC.xlsx'))
sample_data <- read_excel(here('data/sample_annotation/sample_data.xlsx'))


wes_qc <- left_join(wes_qc, sample_data)
wes_qc <- wes_qc %>% filter(!is.na(Type))
passed_wes_qc <- wes_qc %>% filter(WES_STATUS == "PASSED", Type == "Tumor")

g <- ggplot_samplelvl_boxplot(wes_qc, x_variable = "Race", y_variable = "Avg_Cov", 
                              rotate = T, default.theme = "default", legend.pos = "none",
                              y = "Average Coverage", x = "Race")
ggsave(plot = g, path = output.path,  filename = paste0("supfig2Bleft_AvgCov_Race_all.pdf"), device = "pdf", width = 2,  height =3, unit = "in")

g <- ggplot_samplelvl_boxplot(passed_wes_qc, x_variable = "Race", y_variable = "Avg_Cov", 
                              rotate = T, default.theme = "default", legend.pos = "none",
                              y = "Average Coverage", x = "Race")
ggsave(plot = g, path = output.path,  filename = paste0("supfig2Bright_AvgCov_Race_Passed.pdf"), device = "pdf", width = 2,  height =3, unit = "in")


#SupFig1C ####
# Load data
tpm <- readRDS(here("output/01.0_RNAseq_Quality_Control/tpm_all.rds"))
rna_qc <- read_csv(here('output/01.0_RNAseq_Quality_Control/rnaseq_qc.csv'))

keep <- rowSums(tpm) != 0
log2tpm <- log2(tpm[keep,rna_qc$sampleID] + 1)
colMedians(log2tpm) == 0

rna_qc$medlog2tpm0 <- colMedians(log2tpm) == 0
rna_qc <- rna_qc %>% mutate(RNAseq_STATUS = 
                      if_else(medlog2tpm0, "FAILED", "PASSED"))
passed_rna_qc <- rna_qc %>% filter(RNAseq_STATUS == "PASSED")


g <- ggplot_samplelvl_boxplot(passed_rna_qc, x_variable = "Race", y_variable = "total_mapped", 
                              rotate = T, default.theme = "default", legend.pos = "none",y = "Mapped Reads", x = "Race")
ggsave(plot = g, path = output.path,  filename = paste0("PASSED_Race_vs_Mapped_Reads.pdf"), device = "pdf", width = 2,  height =3, unit = "in")

g <- ggplot_samplelvl_boxplot(rna_qc, x_variable = "Race", y_variable = "total_mapped",
                              rotate = T, default.theme = "default", legend.pos = "none",y = "Mapped Reads", x = "Race")
ggsave(plot = g, path = output.path,  filename = paste0("ALL_Race_vs_Mapped_Reads.pdf"), device = "pdf", width = 2,  height =3, unit = "in")

## Plot D (Hierarchical Clustering of Ancestry)
# Load Sample Data. WES and RNAseq QC metrics, and Ancestry ####
sample_data <- read_excel(here('data/sample_annotation/sample_data.xlsx'))
ancestry <- read_csv(here('data/dna/admixture/ancestry_estimates.csv'))

sample_data %>% dplyr::count(ptID) %>% filter(n == 1) #9 patients with Tumor only
sample_data %>% dplyr::count(ptID) %>% filter(n == 2) #57 patients with matched T/N
sample_data %>% select(ptID, Race) %>% unique() %>% dplyr::count(Race)

## WES PIPELINE ####
wes_qc <- wes_qc %>% filter(sampleID %in% sample_data$sampleID)
wes_qc %>% group_by(Group) %>% dplyr::count(WES_STATUS)

wes_paired <- read_excel(here('data/dna/wes_qc/QC.xlsx'), sheet = 2)
wes_tumor_only <- read_excel(here('data/dna/wes_qc/QC.xlsx'), sheet = 3)

wes_paired <- wes_paired %>% filter(`Paired Analysis` == "PASSED")
wes_tumor_only <- wes_tumor_only %>% filter(`Tumor Only` == "PASSED")

wes_qc <- wes_qc %>% mutate(WES_comparison = case_when(
  sampleID %in% wes_paired$sampleID ~ "Matched Normal",
  sampleID %in% wes_tumor_only$sampleID ~ "Pooled Normal",
  .default = NA
))

wes_qc %>% group_by(Group) %>% dplyr::count(WES_comparison)

# SUPFIG 1D ####
plot <- filter(final.cohort, Type == "Tumor")

source(here('code/functions/heatmap_cluster.R'))
source(here('code/functions/heatmap_default.R'))

## Generate Matrix
mat <- as.matrix(plot %>% filter(!is.na(AFR)) %>% select(SAS:AFR))
rownames(mat) <- plot %>% filter(!is.na(AFR)) %>% pull(sampleID)
race <- plot %>% filter(!is.na(AFR))

mat <- t(mat)


## Plot Heatmap
barplot.cols <- c("SAS" = "purple","AMR" = "yellow","EAS" = "orange", "EUR"= "blue","AFR" = "red"  )
library(ComplexHeatmap)
top_anno <- ComplexHeatmap::HeatmapAnnotation(
  height = unit(2, "cm"),
  simple_anno_size = unit(0.2, "cm"),
  show_legend = TRUE,
  annotation_name_gp = gpar(fontsize = 6),
  Ancestry = anno_barplot((t(mat)), gp = gpar(fill = barplot.cols, col = barplot.cols)),
  Race = race$Race,
  col = list(
    Ancestry = barplot.cols,
    Race = c("B" = "red", "W" = "blue")))

dend <- heatmap_cluster(mat)
h <- heatmap_default(mat, coldendogram = dend$col_dend, 
                     title = "% Anc", n.col.branch.col  = 2,
                     top_anno = top_anno, column_split = 2, 
                     show_row_dend = F,  show_column_names = F, 
                     height = unit(0.5, "cm"))

ancestry_legend = Legend(labels = names(barplot.cols), title = "Super Population", legend_gp = gpar(fill = barplot.cols),
                         grid_height = unit(2, "mm"),
                         grid_width = unit(2, "mm"), labels_gp = gpar(fontsize = 4),
                         ncol = 1)


pdf(file = here(output.path, "JHU_Ancestry_clust.pdf"),
    width = 5, height = 3)
ComplexHeatmap::draw(h,  annotation_legend_list = ancestry_legend)
dev.off()

## Cluster
cluster <- cutree(dend$col_dend, k = 2)

## Assign Ancestry Groups
plot <- plot %>% mutate(
  cluster = if_else(
    sampleID %in% names(cluster), cluster[sampleID], NA
  )
)

plot <- plot %>% 
  mutate(ancestry_group = case_when(
    cluster == 1 ~ "AFR", 
    cluster == 2 ~ "EUR"))

## Annotate Final Cohort with Ancestry Group
final.cohort <- left_join(final.cohort, select(plot, ptID, ancestry_group))

# Note, one case without ancestry estimate because sample failed. 

pt.lvl <- final.cohort %>% filter(Type == 'Tumor')
pt.lvl %>% dplyr::count(ancestry_group)
write_csv(pt.lvl, here(output.path, "suptable1.csv"))
write_csv(final.cohort, here(output.path, "sample_lvl.csv"))

