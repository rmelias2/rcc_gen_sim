#' Function to generate histogram or density plot for a given gene (or cytoband) list. 
#' @param diff_genes DESEQ2 Differentially expressed Genes output
#' @param output.path Folder in which to save outputs 
#' @param name File name to save
#' @param volcano_adj.pval summarized experiment object 
#' @param volcano_abs_log2fold Significance cuttoffs 


dea_volcano <- function(diff_genes, output.path, name, volcano_adj.pval = 0.05, volcano_abs_log2fold = 1) {
  required_packages <- c("ggrepel")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("Package", pkg, "is required but not installed. Please install it using install.packages()."))
    } else {
      library(pkg, character.only = TRUE)
    }
  }
  
  results_df <- diff_genes %>% mutate(log10_pvalue = -log10(padj))
  label_subset <- results_df %>% filter(abs(log2FoldChange) >= 1, padj <= volcano_adj.pval)
  
  g <- ggplot(results_df, aes(x = log2FoldChange, y = log10_pvalue)) + 
    geom_point(data = subset(results_df, log2FoldChange >= volcano_abs_log2fold & padj < volcano_adj.pval), color = "red", alpha = 0.7, aes(text = Gene_symbol)) + 
    geom_point(data = subset(results_df, log2FoldChange <= -volcano_abs_log2fold & padj < volcano_adj.pval), color = "blue", alpha = 0.7, aes(text = Gene_symbol)) +
    geom_hline(yintercept = -log10(.05), linetype = 3) + geom_vline(xintercept = c(-volcano_abs_log2fold,volcano_abs_log2fold), linetype = 3) + 
    geom_text_repel(data = label_subset, max.overlaps = 30, force = 0.7, aes(x = log2FoldChange, y = log10_pvalue,  label = Gene_symbol), size = 2.5, min.segment.length	= 0.05
    ) + 
    geom_point(data = subset(results_df, abs(padj) >= volcano_adj.pval | abs(log2FoldChange) < volcano_abs_log2fold & padj <volcano_adj.pval), alpha = 0.05) + theme_classic() + 
    #coord_cartesian(xlim = c(-5, 5)) + 
    ylab(expression("-Log"[10]~Adjusted~italic(P))) + xlab(expression(Log[2])) 
  
  ggsave(path = here(output.path), filename = paste0("volcano_",name, ".png"), plot = g, device = "png", width = 4, height = 2.5, units = "in")
}


