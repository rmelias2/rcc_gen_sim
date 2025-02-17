library(here)
library(tidyverse)
library(ggpubr)

here::i_am("code/fig4.R")
output.path <- here("output/fig4")
if (!dir.exists(here(output.path))) dir.create(here(output.path))


## Prepare Data for Xcell Analysis
counts <- readRDS(here("output/dataprep_pooled/bc_log2tpm_pooled_se.rds"))
meta <- as.data.frame(colData(counts))
tpm <- assays(counts)[["bc.log2.tpm"]]
tpm <- 2^tpm - 1
tpm <- as.data.frame(tpm)
tpm$gene <- rownames(tpm)
tpm <- tpm %>% select(gene, everything())
write_csv(tpm, here(output.path, "xcell_input.csv"))
write_tsv(tpm, here(output.path, "xcell_input.tsv"))


## Xcell was performed using: https://comphealth.ucsf.edu/app/xcell
# Plot XCELL outputs 
XCELL <- read_tsv(here("output/fig4/xCell_results.txt")) %>% dplyr::rename(cell_type = `...1`)
colnames(XCELL) <- str_replace_all(colnames(XCELL), pattern = "\\.", replacement = "-")
XCELL.df <- as.data.frame(t(XCELL[,-1]))
colnames(XCELL.df) <- XCELL$cell_type
XCELL.df$ptID <- rownames(XCELL.df)
XCELL.df <- left_join(XCELL.df, meta)

symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", ""))



##Xcell Heatmap (Fig 4A) ####
source(here('code/functions/heatmap_default.R'))
library(readxl)

##Format as matrix

colnames(XCELL) <- str_replace_all(colnames(XCELL), pattern = "\\.", replacement = "-")
mat <- as.matrix(XCELL[,-1])
rownames(mat) <- XCELL$cell_type



## Subset for cells of interest
cells_of_interest <- read_excel(here('output/fig4/xcell_celltypes.xlsx'))
of_interest <- cells_of_interest %>% filter(!is.na(keep)) %>% pull(xcell)
font_face <- cells_of_interest %>% filter(!is.na(keep)) %>% pull(fontface)
cibersort_mat <- mat[of_interest,as.character(meta$ptID)]
cibersort_mat <- t(scale(t(cibersort_mat)))


## create Dataframe for p values: 
XCELL.df <- as.data.frame(t(XCELL[,-1]))
colnames(XCELL.df) <- XCELL$cell_type
XCELL.df$ptID <- rownames(XCELL.df)
XCELL.df <- left_join(XCELL.df, meta)

## Calculate P values: 
results <- data.frame(cell_type = of_interest, kw.p.val = NA)

for(i in 1:length(of_interest)){
  print(of_interest[i])
  formula <- as.formula(paste("`", of_interest[i], "` ~ ancestry_group", sep = ""))
  kw <- XCELL.df %>% rstatix::kruskal_test(formula = formula)
  results$kw.p.val[i] <- kw$p
}

results$p.adj <- p.adjust(results$kw.p.val, method = "fdr")
results <- results %>% mutate(fdr = case_when(
  p.adj < 0.0001 ~ "****",
  p.adj > 0.0001 & p.adj < 0.001 ~ " ***",
  p.adj > 0.001 & p.adj < 0.01 ~ "  **",
  p.adj > 0.01 & p.adj < 0.05 ~ "   *",
  p.adj > 0.05  ~ " ",
))

library(ComplexHeatmap)
pvals = rowAnnotation(pval = anno_text(results$fdr, gp = gpar(fontsize = 4, fontface = "bold")))




#color_scale = c(-1,1),
meta$sampleID <- meta$ptID
p <- heatmap_default(input.matrix = cibersort_mat, input.meta = meta, title = "Scaled ES",  column_split = meta$ancestry_group,
                     show_column_names = F, show_row_names = T, show_row_dend = T, show_column_dend = F,
                     cluster_column_slices = FALSE, column_title_gp = gpar(fontsize = 8), column_title_rot = 30,
                     left_annotation = pvals, color_scale = c(-1,1),
                     row_names_gp = gpar(fontsize = 4, fontface = font_face))

pdf(file = here(output.path, "XCELL_IMMUNE_heatmap.pdf"), width = 4, height = 4)
ComplexHeatmap::draw(p)
dev.off()


## Fig 4B, 4C Cell Types of Interest

of_interest <- XCELL$cell_type
p = XCELL.df %>% pivot_longer(cols = all_of(of_interest), names_to = "Cell_type", values_to = "Prop")
sig_cells <- results %>% filter(p.adj <= 0.05) %>% pull(cell_type)

cd4t_cells <- 
  c("CD4+ T-cells", "CD4+ memory T-cells",  "Th1 cells", "Th2 cells", 
    "CD4+ naive T-cells", "Tgd cells")
cd8t_cells <- c("CD8+ T-cells", "CD8+ Tcm", "CD8+ Tem", "CD8+ naive T-cells")

dendridic_cell <- c("DC", "aDC", "iDC", "cDC", "pDC")
macrophages <- c("Macrophages", "Macrophages M1" , "Macrophages M2")
b_cell <- c("B-cells", "pro B-cells","naive B-cells" ,"Memory B-cells" , "Class-switched memory B-cells","Plasma cells")

plot_cell_subset <- function(plot){
  ggplot(plot, aes(x = ancestry_group, y = Prop)) + 
    geom_boxplot(aes(fill = ancestry_group), outlier.alpha = 0) + scale_fill_brewer(palette = "Set1") + geom_jitter(size = 0.25, width = 0.1) + 
    facet_wrap(~Cell_type, strip.position = "top", scales = "free", nrow =1) +
    theme_bw() + labs(x = "Ancestry Group", y = "xCELL Enrichment Score") + 
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), 
          text = element_text(size = 6)) +
    guides(fill = FALSE) +
    stat_compare_means(size = 3,
                       label.y.npc = "top",  label = "p.signif", vjust = 0.5, symnum.args = symnum.args)
  
}


p.subset <- p %>% filter(Cell_type %in% sig_cells)
p.subset$Cell_type <- factor(p.subset$Cell_type, levels = c(sig_cells))
g <- plot_cell_subset(p.subset)
ggsave(here(output.path, paste0("X_CELL_signif", "_Tumor_only_AFRvsEUR.pdf")), plot = g, device = "pdf", width = 10, height = 1.5, units = "in")




