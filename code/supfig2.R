library(here)
library(tidyverse)
library(SummarizedExperiment)
library(readxl)
library(MatchIt)
library(Gmisc)

here::i_am("code/supfig2.R")
output.path <- here("output/supfig2")
if (!dir.exists(here(output.path))) dir.create(here(output.path))

# Load GDC KIRC Clinical Manifest ####
clinical <- read_rds(here('data/tcga/clinical_bcr_biotabs.rds'))
names(clinical)
meta <- clinical[["clinical_patient_kirc"]]

## Reformat
meta <- meta[-c(1:2),]
meta <- meta %>% dplyr::rename(ptID = bcr_patient_barcode)


# Identify patients with WXS, RNAseq, CNA, and Ancestry Calls ####
maf <- read_rds(here('data/tcga/kirc_maf.rds'))
maf <- maf %>%
  mutate(
    ptID = sub("^(([^-]*-){2}[^-]*).*", "\\1", Tumor_Sample_Barcode),
    sampleID = sub("^(([^-]*-){3}[^-]*).*", "\\1", Tumor_Sample_Barcode),
    analyteID = str_sub(Tumor_Sample_Barcode, 20, 20)
  )
maf %>% dplyr::count(analyteID)

cns <- read_rds(here('data/tcga/KIRC_DNAcopy_segment_CN.rds'))
cns <- cns %>% dplyr::rename(caseID = Sample)
cns <- cns %>%
  mutate(
    ptID = sub("^(([^-]*-){2}[^-]*).*", "\\1", caseID),
    sampleID = sub("^(([^-]*-){3}[^-]*).*", "\\1", caseID)
  )

se.rna <- read_rds(here('data/tcga/kirc_rnaseq_se.rds'))
rna.meta <- as.data.frame(colData(se.rna))
rna.meta <- rna.meta %>% dplyr::rename(ptID = patient, caseID = barcode, sampleID = sample)

intersect <- BiocGenerics::intersect(maf$ptID, cns$ptID)
intersect <- BiocGenerics::intersect(intersect, rna.meta$ptID)

# Ancestry Information 

#' Ancestry estimations were downloaded from the TCGA publication [Carrot-Zhang et al., Cancer Cell, 2020](https://pubmed.ncbi.nlm.nih.gov/32396860/)
meta.ancestry <- read_xlsx(here("data/tcga/tcga_ancestry.xlsx"), skip = 1) %>% dplyr::rename(ptID = patient)
meta.ancestry %>% dplyr::count(consensus_ancestry)

intersect <- BiocGenerics::intersect(intersect, meta.ancestry$ptID)

#' Note: 356 patients with WXS, RNAseq, CN (DNAcopy_segment), and Ancestry Calls

# Subset meta for patients with all annotation types
meta <- meta %>% filter(ptID %in% intersect)

# Join Ancestry Calls to clinical data 
meta <- left_join(meta, meta.ancestry, by = "ptID")


# Pull Histological Annotation ####

# KIPAN Study 
#' [Ricketts et al., Cell Reports, 2018](https://www.sciencedirect.com/science/article/pii/S2211124718304364?via%3Dihub)
meta.kipan <- read_excel(here("data/tcga/metadata_KIPAN_study.xlsx"))
meta.kipan <- meta.kipan %>% dplyr::rename(ptID = bcr_patient_barcode, 
                                           pankidney_path = `PanKidney Pathology`,
                                           pankidney_pathologic_stage = pathologic_stage,
                                           pankidney_pathologic_T = pathologic_T)


meta <- left_join(meta, dplyr::select(meta.kipan, ptID, pankidney_path, pankidney_pathologic_stage,pankidney_pathologic_T))
meta <- meta %>% mutate(pankidney_path = if_else(is.na(pankidney_path), "missing", pankidney_path))
# Bakouny Study 
#' tRCC cases were identified from supplemental Table 1 of [Bakouny et al., Cell Reports, 2022](https://www.sciencedirect.com/science/article/pii/S2211124721016910?via%3Dihub#app2)
meta.tRCC <- read_xlsx(here("data/tcga/metadata_tRCC_study.xlsx")) %>% dplyr::rename(ptID = `Patient ID`, bak_histology = Histology)
meta.tRCC <- meta.tRCC %>% filter(bak_histology != "Normal") 

meta <- left_join(meta, select(meta.tRCC, ptID, bak_histology))


# Filter nccRCC cases ####
meta %>% dplyr::count(bak_histology, race)
meta %>% dplyr::count(pankidney_path, race)

meta <- meta %>% filter(bak_histology != "tRCC", pankidney_path != "ChRCC")


## 341 Cases remain

# Prepare Propensity Matched Subset ####

#Subset for AFR and EUR ancestries
meta %>% dplyr::count(consensus_ancestry)
meta %>% dplyr::count(race)

meta <- meta %>% filter(consensus_ancestry %in% c("afr", "afr_admix", "eur", "eur_admix"))

meta <- meta %>% dplyr::rename(
  AFR = `Admixture % AFR`,
  EUR = `Admixture % EUR`,
  AMR = `Admixture % AMR`,
  SAS = `Admixture % SAS`,
  EAS = `Admixture % EAS`
)

meta <- meta %>% mutate(predom_ancestry = if_else(consensus_ancestry %in% c("afr", "afr_admix"), "AFR", "EUR"))

# Create Variables on which to match 

## Grade
meta <- meta %>% mutate(tumor_grade = factor(if_else(tumor_grade == "[Not Available]", "GX",tumor_grade)))

## Stage
meta <- meta %>% mutate(stage = factor(if_else(ajcc_pathologic_tumor_stage == "[Discrepancy]", "NA",ajcc_pathologic_tumor_stage)))

## Gender 
meta$gender <- factor(meta$gender)

## Ancestry
meta$predom_ancestry <- factor(meta$predom_ancestry, levels = c("EUR", "AFR"))


### Matching ####
m.out <- matchit(predom_ancestry ~ gender + ajcc_pathologic_tumor_stage + tumor_grade, data = meta, method = "optimal", distance = "glm", 
                 ratio = 3)
summary(m.out)

matched.kirc.meta <- match.data(m.out)

### Visualize Propensity Matched Output ####
matched.kirc.meta <- matched.kirc.meta %>% mutate(age_at_initial_pathologic_diagnosis = as.numeric(age_at_initial_pathologic_diagnosis))

matched.kirc.meta <- matched.kirc.meta %>% mutate(race = case_when(
  race == "BLACK OR AFRICAN AMERICAN" ~ "B",
  race == "WHITE" ~ "W",
  .default = "NA"
),
race = factor(race, levels = c("B", "W", "NA")))


matched.kirc.meta <- matched.kirc.meta %>% mutate(gender = factor(gender)) %>% 
  set_column_labels(
    age_at_initial_pathologic_diagnosis = "Age",
    gender = "Sex",
    tumor_grade = "Nuclear Grade",
    stage = "Stage",
    race = "Race",
    predom_ancestry = "Consensus Ancestry"
  ) 


t <- matched.kirc.meta %>% 
  getDescriptionStatsBy(age_at_initial_pathologic_diagnosis, gender, tumor_grade, stage	, race,
                        by = predom_ancestry, header_count = TRUE, show_all_values = TRUE, statistics = T)  %>% 
  htmlTable()

htmltools::save_html(htmlTable(t), file = here(output.path, "tcga_matched_cohort.html"))


# Prepare Metadata File for Downstream Analyses ####

matched.kirc.meta.clean <- matched.kirc.meta %>% dplyr::select(ptID, race, gender, age_at_initial_pathologic_diagnosis, stage, ajcc_tumor_pathologic_pt, tumor_grade, predom_ancestry, 
                                                               consensus_ancestry, AFR:SAS)
matched.kirc.meta.clean <- matched.kirc.meta.clean %>% dplyr::rename(Sex = gender, Stage = stage, pathologic_t = ajcc_tumor_pathologic_pt,
                                                                     path_grade = tumor_grade, 
                                                                     Race = race)

colnames(matched.kirc.meta.clean)


matched.kirc.meta.clean <- matched.kirc.meta.clean %>% mutate(
  pathologic_t_group = case_when(
    str_detect(pathologic_t, "T1") ~ "pT1",
    str_detect(pathologic_t, "T2") ~ "pT2",
    str_detect(pathologic_t, "T3") ~ "pT3",
    str_detect(pathologic_t, "T4") ~ "pT4"
  ),
  path_grade = str_replace(path_grade, "G", "")
)


# Assign Ancestry Group ####
# Ancestry group is assigned via clustering of % ancestries 

source(here('code/functions/heatmap_cluster.R'))
source(here('code/functions/heatmap_default.R'))

## Generate Matrix
mat <- as.matrix(matched.kirc.meta.clean %>% filter(!is.na(AFR)) %>% select(AFR:SAS))
rownames(mat) <- matched.kirc.meta.clean %>% filter(!is.na(AFR)) %>% pull(ptID)

race <-  matched.kirc.meta.clean %>% filter(!is.na(AFR))

mat <- t(mat)


## Plot Heatmap
barplot.cols <- c("AFR" = "red","AMR" = "yellow", "EAS" = "orange",  "EUR"= "blue", "SAS" = "purple")
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
    Race = c("B" = "red", "W" = "blue", "NA" = "white"))
)

dend <- heatmap_cluster(mat)
h <- heatmap_default(mat, coldendogram = dend$col_dend, title = "% Anc",  
                     top_anno = top_anno, column_split = 2, show_row_dend = F,  
                     show_column_names = F, height = unit(1, "cm"))

ancestry_legend = Legend(labels = names(barplot.cols),
                         title = "Super Population", 
                         legend_gp = gpar(fill = barplot.cols),
                         grid_height = unit(2, "mm"),
                         grid_width = unit(2, "mm"), labels_gp = gpar(fontsize = 4),
                         ncol = 1)


pdf(file = here(output.path, "TCGA_Ancestry_clust.pdf"),
    width = 5, height = 3)
ComplexHeatmap::draw(h,  annotation_legend_list = ancestry_legend)
dev.off()

## Cluster
cluster <- cutree(dend$col_dend, k = 2)

## Assign Ancestry Groups
matched.kirc.meta.clean <- matched.kirc.meta.clean %>% mutate(
  cluster = if_else(
    ptID %in% names(cluster), cluster[ptID], NA
  )
)


matched.kirc.meta.clean <- matched.kirc.meta.clean %>% 
  mutate(ancestry_group = case_when(
    cluster == 1 ~ "AFR", 
    cluster == 2 ~ "EUR")) %>% 
  mutate(ancestry_group = if_else(is.na(ancestry_group), as.character(predom_ancestry), ancestry_group))

# Save Metadata File ####
matched.kirc.meta.clean %>% write_tsv(here(output.path, "matched_kirc_meta.tsv"))


