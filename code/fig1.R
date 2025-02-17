library(here)
library(tidyverse)
library(ComplexHeatmap)

here::i_am("code/fig1.R")
output.path <- here("output/fig1")
if (!dir.exists(here(output.path))) dir.create(here(output.path))

source(here('code/functions/oncoprint_default.R'))
source(here('code/functions/oncoprint_split.R'))
source(here('code/functions/format_maf_to_matrix.R'))

## load data

maf <- read_tsv(here('output/dataprep_jhu/jhu_cohort.maf'))
## Note, due to patient privacy concerns, MAF files are unable to be uploaded with this manuscript. 
## Please provide the corresponding author a request for data and we will share if appropriate safety measures are in place.
meta <- read_csv(here('output/supfig1/suptable1.csv'))
samples <- meta %>% filter(WES_STATUS != "FAILED") %>% pull(sampleID)


## Prepare gene matrix
goi <- c("VHL", "BAP1", "SETD2", "PBRM1",
         "PTEN", "TSC1", "TSC2", "CDKN2A","KDM5C", "MTOR","TP53","HIF1A")

mat <- format_maf_to_matrix(maf.df = maf, genes = goi)

mat <- annotate_no_call(mat, samples)


meta <- meta %>% filter(sampleID %in% colnames(mat))
mut.mat <- mat[,meta$sampleID]


meta$ancestry_group <- factor(meta$ancestry_group, levels = c("EUR", "AFR"))

## Figure 1A Driver Mutations in ccRCC by ancestry group
p <- oncoprint_split(input.matrix = mut.mat, 
                     input.meta = meta, 
                     split_var = "ancestry_group", 
                     height = unit(1.3, "in"),
                     show_column_names = FALSE)

pdf(file = here(output.path, "goi_oncoprint.pdf"),
    width = 5, height = 2)
ComplexHeatmap::draw(p)
dev.off()


## Figure 1B, C, D

source(here('code/functions/format_maf_to_matrix.R'))
source(here('code/functions/stats_group_frequencies_and_pvals.R'))
source(here('code/functions/ggplot_alteration_frequency_scatter.R'))

maf <- read_tsv(here('output/dataprep_jhu/jhu_cohort.maf'))

vhl.mut <- maf %>% filter(Hugo_Symbol == "VHL") 

meta <- meta %>% mutate(vhl = if_else(sampleID %in% vhl.mut$Tumor_Sample_Barcode, "Mut", "Wt"))

meta$vhl <- factor(meta$vhl, levels = c("Wt", "Mut"))

source(here('code/functions/ggplot_stackedbar.R'))
cols = c("Mut" = "darkred", "Wt" = "gray")


g <- ggplot_stackedbar(meta, x_variable = "Race", group_variable = "vhl", group.colors = cols, fill = "VHL Status", x = "Race", y = "Frequency", annotate_n = T)
ggsave(plot = g, filename = "vhl_Race.pdf", path = here(output.path), width = 3, height = 4, device = "pdf")

g <- ggplot_stackedbar(meta, x_variable = "HTN", group_variable = "vhl", fill = "VHL Status", group.colors = cols, x = "HTN", y = "Frequency", annotate_n = T)
ggsave(plot = g, filename = "vhl_HTN.pdf", path = here(output.path), width = 3, height = 4, device = "pdf")

g <- ggplot_stackedbar(meta, x_variable = "GFR_bin", group_variable = "vhl", fill = "VHL Status", x = "GFR", group.colors = cols, y = "Frequency", annotate_n = T)
ggsave(plot = g, filename = "vhl_GFR.pdf", path = here(output.path), width = 3, height = 4, device = "pdf")

