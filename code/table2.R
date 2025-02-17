library(here)
library(SummarizedExperiment)
library(tidyverse)

here::i_am("code/table2.R")
output.path <- here("output/table2")
if (!dir.exists(here(output.path))) dir.create(here(output.path))


## Loading JHU Cohort
jhu.meta <- read_csv(here('output/supfig1/suptable1.csv'))
jhu.meta$cohort <- "JHU"
jhu.meta <- jhu.meta %>% filter(Type == "Tumor")

tcga.meta <- read_tsv(here('output/supfig2/matched_kirc_meta.tsv'))
tcga.meta$cohort <- "TCGA"
tcga.meta$RNAseq_STATUS <- "PASSED"
tcga.meta$WES_comparison <- "Pooled Normal" 
tcga.meta$WES_STATUS <- "PASSED"
colnames(jhu.meta)
colnames(tcga.meta)
tcga.meta <- tcga.meta %>% dplyr::rename(Age = age_at_initial_pathologic_diagnosis)

overlap <- intersect(colnames(jhu.meta), colnames(tcga.meta))
meta <- rbind(select(jhu.meta, all_of(overlap)), select(tcga.meta, all_of(overlap)))


library(Gmisc)
meta <- meta %>% mutate(Stage = case_when(
  Stage == "Locally Advanced" ~ "Stage III",
  Stage == "pT1a" ~ "Stage I",
  .default = Stage)
)
meta <- meta %>% mutate(Sex = case_when(
  Sex == "FEMALE" ~ "Female",
  Sex == "MALE" ~ "Male",
  .default = Sex
))


# Patient Level

df <- meta
df$path_grade <- factor(df$path_grade)
df$Sex  <- factor(df$Sex)
df$Stage <- factor(df$Stage)

df <- df %>% set_column_labels(
  path_grade = "Grade",
  ancestry_group = "Estimated Ancestry",
  cohort = "Cohort",
  WES_STATUS = "WES Available",
  RNAseq_STATUS = "RNAseq Available"
) 


t <- df %>% getDescriptionStatsBy(Age, Sex, Stage, path_grade, cohort, WES_STATUS, RNAseq_STATUS,
                                  by = ancestry_group, add_total_col = TRUE, NEJMstyle = T,
                                  header_count = TRUE, show_all_values = T, statistics = T) %>% 
  htmlTable()

htmltools::save_html(htmlTable(t), file = here(output.path, paste0("pooled_cohort_baseline", ".html")))

meta %>% write_csv(here(output.path, "pooled_cohort.csv"))



