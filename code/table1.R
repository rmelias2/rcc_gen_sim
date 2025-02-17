library(here)
library(tidyverse)
library(Gmisc)


here::i_am("code/table1.R")
output.path <- here("output/table1")
if (!dir.exists(here(output.path))) dir.create(here(output.path))

# Table 1 ####
data <- read_csv(here('output/supfig1/suptable1.csv'))
pt.lvl <- data %>% filter(Type == "Tumor")
pt.lvl <- pt.lvl %>% mutate(pT = if_else(pathologic_t == "pT1a", "pT1", "pT3"))
pt.lvl <- pt.lvl %>% mutate(WES_comparison = if_else(is.na(WES_comparison), "Failed",WES_comparison))
pt.lvl <- pt.lvl %>% mutate(RNAseq_STATUS = if_else(RNAseq_STATUS == "FAILED", "Failed", "Passed"))

# Patient Level
pt.lvl$path_grade <- factor(pt.lvl$path_grade)
pt.lvl$WES_comparison <- factor(pt.lvl$WES_comparison, levels = c("Matched Normal", "Pooled Normal", "Failed"))

df <- pt.lvl
df$Sex  <- factor(df$Sex)
df$pT <- factor(df$pT)
df$HTN <- factor(df$HTN)
df$GFR_bin <- factor(df$GFR_bin)
df$RNAseq_STATUS <- factor(df$RNAseq_STATUS, levels = c("Passed", "Failed"))

df <- df %>% set_column_labels(
  pT = "pT Stage",
  year_group = "Year of Surgery",
  WES_comparison = "WES Available",
  path_grade = "Grade",
  RNAseq_STATUS = "RNAseq Available",
  ancestry_group = "Estimated Ancestry",
  eGFRcr = "eGFR",
  GFR_bin = "GFR Group"
) 

df %>% write_csv(here(output.path, "suptab1.csv"))

t <- df %>% getDescriptionStatsBy(Age, Sex, eGFRcr, HTN, pT, path_grade, ancestry_group, RNAseq_STATUS, WES_comparison, 
                                  by = Race, add_total_col = TRUE, NEJMstyle = T,
                                  header_count = TRUE, show_all_values = T, statistics = T) %>% 
  htmlTable()

htmltools::save_html(htmlTable(t), file = here(output.path, paste0("ptlvl", ".html")))
