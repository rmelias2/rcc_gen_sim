library(here)
library(tidyverse)
library(readxl)
library(ComplexHeatmap)

here::i_am("code/fig5.R")
output.path <- here("output/fig5")
if (!dir.exists(here(output.path))) dir.create(here(output.path))


## Load Data
se <- read_rds(here('output/dataprep_pooled/bc_log2tpm_pooled_se.rds'))
meta <- as.data.frame(colData(se))

## Assign Color Parameters for Cluster Subtypes
IMM151_Predictions <- read_csv(here('output/supfig4/IMM151_predications.csv'))
meta <- left_join(meta, select(IMM151_Predictions, ptID, cluster))

meta <- meta %>% 
  mutate(cluster_label = case_when(
    cluster == "angiogenic" ~ "Angiogenic",
    cluster ==   "angio_stromal" ~ "Angio/Stromal",
    cluster ==   "proliferative" ~ "Proliferative",
    cluster ==   "complement_omega_oxidation" ~ "Complement/Omega-ox.",
    cluster ==   "teff_proliferative" ~ "T-eff/Proliferative",
    cluster ==   "stromal_proliferative" ~ "Stromal/Proliferative")
  )

meta$cluster_label <-
  factor(
    meta$cluster_label,
    levels = c(
      "Angiogenic",
      "Angio/Stromal",
      "Complement/Omega-ox.",
      "T-eff/Proliferative",
      "Proliferative",
      "Stromal/Proliferative"
    )
  )

subtype.col <- c("Angiogenic" = "#6d1f82",
                 "Angio/Stromal" =  "#fa170f",
                 "Stromal/Proliferative" = "gray3",
                 "Proliferative" = "#faa80f", 
                 "T-eff/Proliferative" = "#1f8248", 
                 "Complement/Omega-ox." = "#539bed")

meta %>% write_rds(here(output.path, "IMM151_clusters.rds"))
## Pull Genes of Interest
IMm151_genes <- read_excel(here('data/gene_signatures/rcc_gene_signatures.xlsx'))
IMm151_genes <- IMm151_genes %>% filter(str_detect(Source, "Motzer"), Label != "snoRNA")

IMm151_genes$Label <- factor(
  IMm151_genes$Label,
  levels = c(
    "T-effector",
    "Angiogenesis",
    "FAO AMPK",
    "Cell Cycle",
    "FAS Pentose Phosphate",
    "Stroma",
    "Myeloid Inflammation",
    "Complement Cascade",
    "Omega-Oxidation"
  )
)

source(here('code/functions/heatmap_default.R'))
source(here('code/functions/heatmap_top_annotation_colors.R'))

mat <- assays(se)[["bc.log2.tpm"]]
mat <- mat[IMm151_genes$Genes,]
mat <- t(scale(t(mat)))

top_annotation = HeatmapAnnotation(
  annotation_legend_param = annotation_legend_param,
  show_legend = TRUE,
  simple_anno_size = unit(0.2, "cm"),
  annotation_name_gp = gpar(fontsize = 5.5),
  AFR = meta$AFR,
  EUR = meta$EUR,
  Subtype = meta$cluster_label,
  col = list(AFR = AFR.col, EUR = EUR.col,
             Subtype = subtype.col))

meta$sampleID <- meta$ptID
h <- heatmap_default(input.matrix = mat, input.meta = meta, title = "Normalized Exp.", 
                     show_column_dend = F, show_column_names = F, split = IMm151_genes$Label,
                     column_split = meta$cluster, row_title_rot = 0, column_title = " ",
                     cluster_column_slices = F,
                     top_anno = top_annotation, row_title_gp = gpar(fontsize = 6, fontface = "bold"))

pdf(file = here(output.path, "IM151_Heatmap_clust.pdf"),
    width = 5, height = 5)
ComplexHeatmap::draw(h)
dev.off()


## Stacked Bar by Subgroup


source(here('code/functions/ggplot_stackedbar.R'))

##Subtype and Ancestry 
g <- ggplot_stackedbar(meta, x_variable = "ancestry_group", 
                       group_variable = "cluster_label", 
                       fill = "IMM151 Subtype", x = "Ancestry Group", 
                       y = "Frequency", stat = "fisher.test", group.colors = subtype.col, pval.x.pos = 1, annotate_n = T)
g <- g + theme(legend.position = "none")
ggsave(plot = g, filename = "IMM151_Ancestry_Bin.pdf", path = here(output.path), width = 3, height = 5, device = "pdf")


## GSVA Analysis: 
library(readxl)
library(GSVA)
library(fgsea)


se <- read_rds(here('output/dataprep_pooled/bc_log2tpm_pooled_se.rds'))
meta <- as.data.frame(colData(se))

log_tpm <- assays(se)[["bc.log2.tpm"]]

#Load gene sets
gene_sigs <- read_excel(here("data/gene_signatures/rcc_gene_signatures.xlsx"))
gene_sets <-
  (gene_sigs %>% group_by(Classification) %>% summarise(genes = list(unique(Genes))) %>% as.list(.))[[2]]
names(gene_sets) <-
  (gene_sigs %>% group_by(Classification) %>% summarise(genes = list(unique(Genes))) %>% as.list(.))[[1]]

HALLMARK <- gmtPathways(here("data/gene_signatures/msigdb_v2022.1.Hs_GMTs/h.all.v2022.1.Hs.symbols.gmt"))
gene_sets <- c(gene_sets, HALLMARK)

#GSVA
param.gsva <- gsvaParam(
  exprData = log_tpm,
  geneSets = gene_sets,
  minSize = 4
)


gsva.output<- gsva(param.gsva)

#Output
gsva.df <- as.data.frame(t(gsva.output))
gsva.df$ptID <- colnames(gsva.output)
meta <- as.data.frame(colData(se))
meta <- left_join(meta, gsva.df)

meta %>% write_csv(here(output.path, "gsva_results.csv"))


# Visualize: 


source(here('code/functions/ggplot_samplelvl_boxplot.R'))

# LOAD DATA and SET GRAPHICS ####
subtype.col <- c("Angiogenic" = "#6d1f82",
                 "Angio/Stromal" =  "#fa170f",
                 "Stromal/Proliferative" = "gray3",
                 "Proliferative" = "#faa80f", 
                 "T-eff/Proliferative" = "#1f8248", 
                 "Complement/Omega-ox." = "#539bed")

symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", ""))


gsva <- read_csv(here('output/fig5/gsva_results.csv'))
IMM151 <- read_rds(here('output/fig5/IMM151_clusters.rds'))


plot <- left_join(gsva, select(IMM151, ptID, cluster, cluster_label))


colnames(plot) <- colnames(plot) %>% str_remove("HALLMARK_")



select.genesets <- c("ALLOGRAFT_REJECTION", "BILE_ACID_METABOLISM", "E2F_TARGETS", "FATTY_ACID_METABOLISM","G2M_CHECKPOINT", "MYC_TARGETS_V1", "GLYCOLYSIS", "INTERFERON_GAMMA_RESPONSE",
                     "INTERFERON_ALPHA_RESPONSE", "MTORC1_SIGNALING","OXIDATIVE_PHOSPHORYLATION", "GLYCOLYSIS", "PEROXISOME", "UNFOLDED_PROTEIN_RESPONSE")

p = plot %>% pivot_longer(cols = all_of(select.genesets), names_to = "Gene Set", values_to = "GSVA")

library(ggpubr)

# Sup Figure (GSVA of significant hallmarks)
g <- ggplot(p, aes(x = ancestry_group, y = GSVA)) + 
  geom_boxplot(aes(fill = ancestry_group), outlier.alpha = 0) + 
  geom_jitter(size = 0.25, width = 0.1) + 
  scale_fill_manual(values = c("AFR" = "red", "EUR" = "blue")) + 
  facet_wrap(~`Gene Set`, strip.position = "top", scales = "free", nrow = 3) +
  theme_bw() + labs(y = "GSVA", x = "IMM151 Subtype", fill = "IMM151 Subtype") + 
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 6)) +
  stat_compare_means(size = 3,
                     label.y.npc = "top",  label = "p.signif", vjust = 0.5, symnum.args = symnum.args) 

ggsave(plot = g, filename = "SELECT_Hallmark_Ancestry.pdf", path = here(output.path),
       width = 7, height =6)

select.genesets <- c("ALLOGRAFT_REJECTION", "BILE_ACID_METABOLISM", "E2F_TARGETS", "FATTY_ACID_METABOLISM","G2M_CHECKPOINT", "GLYCOLYSIS", "INTERFERON_GAMMA_RESPONSE",
                     "INTERFERON_ALPHA_RESPONSE", "MTORC1_SIGNALING", "GLYCOLYSIS", "PEROXISOME")

p = plot %>% pivot_longer(cols = all_of(select.genesets), names_to = "Gene Set", values_to = "GSVA")


## Fig 5 C ###
g <- ggplot(p, aes(x = cluster_label, y = GSVA, group = interaction(cluster_label, ancestry_group))) + 
  geom_boxplot(aes(fill = ancestry_group), outlier.alpha = 1, outlier.size = 0.25) +   
  scale_fill_manual(values = c("AFR" = "red", "EUR" = "blue")) + 
  facet_wrap(~`Gene Set`, strip.position = "top", scales = "free", nrow = 2) +
  theme_bw() + labs(y = "GSVA", x = "IMM151 Subtype", fill = "Ancestry Group") + 
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none",
        text = element_text(size = 6)) +
  stat_compare_means(size = 3,
                     label.y.npc = "top",  label = "p.signif", vjust = 0.5, symnum.args = symnum.args) 

ggsave(plot = g, filename = "SELECT_Ancestry_Clust_Interaction_HALLMARK.pdf", path = here(output.path),
       width = 8, height = 5)

## Fig 5 D ###
g <- ggplot(p, aes(x = cluster_label, y = GSVA)) + 
  geom_boxplot(aes(fill = cluster_label), outlier.alpha = 0) + 
  geom_jitter(size = 0.25, width = 0.1) + 
  scale_fill_manual(values = subtype.col) + 
  facet_wrap(~`Gene Set`, strip.position = "top", scales = "free", nrow = 2) +
  theme_bw() + labs(y = "GSVA", x = "IMM151 Subtype", fill = "IMM151 Subtype") + 
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        text = element_text(size = 6)) +
  stat_compare_means(size = 3,
                     label.y.npc = "top",  label = "p.signif", vjust = 0.5, symnum.args = symnum.args) 

ggsave(plot = g, filename = "SELECT_CLUST_HALLMARK.pdf", path = here(output.path),
       width = 7, height =2.5)




