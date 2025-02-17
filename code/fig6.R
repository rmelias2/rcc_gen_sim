library(here)
library(tidyverse)
library(SummarizedExperiment)
here::i_am("code/fig6.R")
output.path <- here("output/fig6.R")
if (!dir.exists(here(output.path))) dir.create(here(output.path))

## A - Oncoprint by IMM151 Clust ####
maf <- read_tsv(here('output/dataprep_pooled/pooled_maf.tsv'))
IMM151 <- read_rds(here('output/fig5/IMM151_clusters.rds'))
meta <- read_csv(here('output/table2/pooled_cohort.csv'))


meta <- left_join(meta, select(IMM151, cluster_label, ptID ))
meta <- meta %>% filter(!is.na(cluster_label))

source(here('code/functions/format_maf_to_matrix.R'))


goi <- c("VHL", "BAP1", "SETD2", "PBRM1",
         "PTEN", "TSC1", "TSC2","KDM5C", "MTOR","TP53","HIF1A")

gene_level_mat <- format_maf_to_matrix(maf.df = maf, genes = goi, sampleID = "ptID")
gene_level_mat <- annotate_no_call(gene_level_mat, samples = meta$ptID)


# Plot Oncoprint
source(here('code/functions/oncoprint_default.R'))
source(here('code/functions/oncoprint_split.R'))
library(ComplexHeatmap)

comparison.df <- meta 
comparison.df$sampleID = meta$ptID
mat <- gene_level_mat
comparison.df <- comparison.df %>% filter(sampleID %in% colnames(mat))
mut.mat <- mat[,comparison.df$sampleID]
goi <- rownames(mut.mat)

comparison.df$cluster_label <-
  factor(
    comparison.df$cluster_label,
    levels = c(
      "Angiogenic",
      "Angio/Stromal",
      "Complement/Omega-ox.",
      "T-eff/Proliferative",
      "Proliferative",
      "Stromal/Proliferative"
    )
  )


p <- oncoprint_split(input.matrix = mut.mat, 
                     input.meta = comparison.df, 
                     split_var = "cluster_label", 
                     height = unit(1.3, "in"),
                     show_column_names = FALSE)

pdf(file = here(output.path, "goi_oncoprint.pdf"),
    width = 9, height = 2)
ComplexHeatmap::draw(p)
dev.off()



## Fig 6B
## VHL Mutation by Subtype ####
maf <- read_tsv(here('output/dataprep_pooled/pooled_maf.tsv'))
IMM151 <- read_rds(here('output/fig5/IMM151_clusters.rds'))
meta <- read_csv(here('output/table2/pooled_cohort.csv'))


meta <- left_join(meta, select(IMM151, cluster_label, ptID ))
meta <- meta %>% filter(!is.na(cluster_label))
vhl.mut <- maf %>% filter(Hugo_Symbol == "VHL")
meta <- meta %>% mutate(VHL = if_else(ptID %in% vhl.mut$ptID, "mut", "wt"))
source(here('code/functions/ggplot_stackedbar.R'))
unique(meta$cluster_label)
subtype <- unique(meta$cluster_label)
meta$VHL <- factor(meta$VHL, levels = c("wt", "mut"))
ver1 <- meta %>% filter(cluster_label == "T-eff/Proliferative")

for(i in subtype) {
  ver1 <- meta %>% filter(cluster_label == i)
  g <- ggplot_stackedbar(df = ver1, x_variable = "ancestry_group", group_variable = "VHL", annotate_n = T, x = "Ancestry Group", y = "Frequency", group.colors = c("mut" = "red", "wt" = "gray"))
ggsave(plot = g, path = here(output.path), filename = (paste0(str_replace_all(i, "[./-]", "_") , "_vhl_by_ancestry.pdf")), device = "pdf", width = 3, height = 3)
}
sessionInfo()

