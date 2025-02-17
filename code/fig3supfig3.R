library(here)
library(tidyverse)
library(edgeR)
library(DESeq2)
library(doParallel)
library(tidyverse)
library(ggpubr)
library(ggrepel)

source(here('code/functions/dea_data_prep.R'))
source(here('code/functions/dea_save_output.R'))
source(here('code/functions/dea_volcano.R'))

here::i_am("code/fig3supfig3v2.R")
output.path <- here("output/fig3supfig3v2")
if (!dir.exists(here(output.path))) dir.create(here(output.path))


# Fig 3, Sup Fig 3 ####
## JHU DEA ####
se <- readRDS(here("output/01.2_metadata_prep/rna_se.rds"))
se$ancestry_bin <- factor(se$ancestry_bin, levels = c("EUR", "AFR"))
se$ancestry_type <- paste0(se$ancestry_bin, "_", se$Type)
se$ancestry_type <- factor(se$ancestry_type, levels = c("EUR_Normal","AFR_Normal", "EUR_Tumor", "AFR_Tumor"))

meta <- as.data.frame(colData(se))
meta %>% dplyr::count(ancestry_type)
counts <- dea_data_prep(se.object = se, 
                        group.var = "ancestry_type")

dds <- DESeqDataSet(counts, design = ~ancestry_type)
dds <-  DESeq(dds,  parallel = TRUE)


#Save Results
resultsNames(dds)[-1]
extra.comparisons <- list(
  "AFR_T_vs_AFR_N" = c("ancestry_type", "AFR_Tumor", "AFR_Normal"),
  "AFR_N_vs_EUR_N" =  c("ancestry_type", "AFR_Normal", "EUR_Normal"),
  "EUR_T_vs_EUR_N" = c("ancestry_type", "EUR_Tumor", "EUR_Normal"),
  "AFR_T_vs_EUR_T" = c("ancestry_type", "AFR_Tumor", "EUR_Tumor")
)

dea_save_output(deseq2.output = dds, exp.name = "Ancestry_Bin_wNormal",
                additional.comparisons = extra.comparisons, save.default = F)


## SupFig3A Ancestry TvsN ####
AFR <- read_tsv(here('output/fig3supfig3v2/Ancestry_Bin_wNormal/Results_AFR_T_vs_AFR_N.txt'))
EUR <- read_tsv(here('output/fig3supfig3v2/Ancestry_Bin_wNormal/Results_EUR_T_vs_EUR_N.txt'))

AFR.df <- dplyr::select(AFR,Gene_symbol, stat,padj) %>% dplyr::rename(AFR = stat,AFR.padj = padj)
EUR.df <- dplyr::select(EUR,Gene_symbol, stat, padj) %>% dplyr::rename(EUR = stat, EUR.padj = padj)

plot <- left_join(AFR.df, EUR.df)

# Add Deseq2 Wald test statistic
plot <- plot %>%
  mutate(significance = case_when(
    AFR.padj < 0.05 | EUR.padj < 0.05 ~ "Significant",
    TRUE ~ "Not Significant"
  ))

# Filter for significant points in both ancestry.padj and type.padj
sig_both <- plot %>% filter(significance != "Not Significant")

#Plot
p <- ggplot(plot, aes(x = AFR, y = EUR)) +
  geom_point(size = 1.5, aes(color = significance, alpha = significance)) +  
  scale_color_manual(values = c("Not Significant" = "gray", 
                                "Significant" = "red")) +
  #geom_smooth(method = "lm", color = "black", se = TRUE) +  # Best fit line with confidence interval
  scale_alpha_manual(values = c("Not Significant" = 0.7, 
                                "Significant" = 0.7), guide = "none") +
                     
  annotate("text", x = 12, y = -0.2, label = "AFR T", size = 5, vjust = 1) +
  annotate("text", x = -12, y = -0.2, label = "AFR N", size = 5, vjust = 1) +
  annotate("text", x = 0.5, y = 10, label = "EUR T", size = 5, hjust = 0) +
  annotate("text", x = 0.5, y = -10, label = "EUR N", size = 5, hjust = 0) +
  stat_cor(label.x = -15, label.y = 10) +
  # Add Intercepts
  # geom_vline(xintercept = 0, linetype = "solid", color = "black", size = 0.8) +
  # geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.8) +
  annotate("segment", x = -15, xend = 15, y = 0, yend = 0, color = "black", size = 0.8,
           arrow = arrow(type = "closed", ends = "both", length = unit(0.2, "cm"))) +
  annotate("segment", x = 0, xend = 0, y = -15, yend = 15, color = "black", size = 0.8,
           arrow = arrow(type = "closed", ends = "both", length = unit(0.2, "cm"))) +
  scale_x_continuous(breaks = seq(-15, 15, by = 3)) +
  scale_y_continuous(breaks = seq(-15, 15, by = 3)) +
  
  # Label points
  #geom_text_repel(data = sig_both, aes(label = Gene_symbol), size = 3, max.overlaps = 20, force = 2, guide = "none") +
  labs(x = "Wald Statistic", y = "Wald Statistic", color = "Significance") + 
  # Add axis ticks
  theme_minimal(base_size = 12) + 
  theme(
    axis.line = element_blank(),  
    axis.text = element_text(size = 12),
    panel.grid = element_blank(),   
    panel.border = element_blank()
  )

ggsave(plot = p, path = here(output.path), filename = "wald_AFRTvsN_EURTvsN.png", width = 6, height = 5, dpi = 600, bg = "white")

## Perform TCGA DEA ####

## Load SE
se <- readRDS(here("output/dataprep_TCGA/gene_expression_tcga_se.rds"))
se <- se[,se$Type == "Tumor"]

## Perform DEA
se$ancestry_group <- factor(se$ancestry_group, levels = c("EUR", "AFR"))

meta <- as.data.frame(colData(se))
meta %>% dplyr::count(ancestry_group)

counts <- dea_data_prep(se.object = se, 
                        group.var = "ancestry_group")

dds <- DESeqDataSet(counts, design = ~ancestry_group)
dds <-  DESeq(dds,  parallel = TRUE)

## Save Results
resultsNames(dds)[-1]


dea_save_output(deseq2.output = dds, exp.name = "TCGA_Ancestry_Tumor",
                save.default = T)

tcga <- read_tsv(here('output/fig3supfig3v2/TCGA_Ancestry_Tumor/Results_ancestry_group_AFR_vs_EUR.txt'))


## Identify overlapping AFR T and EUR T associated genes ####


### Load JHU DEA results 
Tumor <- read_tsv(here('output/fig3supfig3/Ancestry_Type_Merged/Results_AFR_T_vs_EUR_T.txt')) # (comparing AFR T vs EUR T )
Normal <- read_tsv(here('output/fig3supfig3/Ancestry_Type_Merged/Results_AFR_N_vs_EUR_N.txt')) # (comparing AFR N vs EUR N)


jhu.genes.up <- Tumor %>% filter(padj < 0.05, stat > 0)
jhu.genes.down <- Tumor %>% filter(padj < 0.05, stat < 0)
tcga.genes.up <- tcga %>% filter(padj < 0.05, stat > 0)
tcga.genes.down <- tcga %>% filter(padj < 0.05, stat < 0)

shared_up <- jhu.genes.up$Gene_symbol[jhu.genes.up$Gene_symbol %in% tcga.genes.up$Gene_symbol]
shared_down <- jhu.genes.down$Gene_symbol[jhu.genes.down$Gene_symbol %in% tcga.genes.down$Gene_symbol]
shared <- c(shared_up, shared_down)



## Fig 3A Plot Interaction JHU cohort ####
Tumor.df <- dplyr::select(Tumor,Gene_symbol, stat,padj) %>% dplyr::rename(Tumor = stat,Tumor.padj = padj)
Normal.df <- dplyr::select(Normal,Gene_symbol, stat, padj) %>% dplyr::rename(Normal = stat, Normal.padj = padj)

plot <- left_join(Tumor.df, Normal.df)

remove <- c("MIR6723", "SNORD115-6") # This has been dropped by NCI gene cards
plot <- plot %>% filter(!Gene_symbol %in% remove)

# Pull Wald Statistic
plot <- plot %>%
  mutate(significance = case_when(
    Tumor.padj < 0.05 & Normal.padj < 0.05 ~ "Both",
    Tumor.padj < 0.05 & Normal.padj >= 0.05 ~ "Tumor Only",
    Tumor.padj >= 0.05 & Normal.padj < 0.05 ~ "Normal Only",
    TRUE ~ "Not Significant"
  ))

# Label genes which are significant in both TCGA and JHU comparisons
to_lab <- shared
sig_both <- plot %>% filter(Gene_symbol %in% to_lab)

DESeq2::plotCounts(dds, "PADI1", intgroup = "ancestry_type")

#Plot
p <- ggplot(plot, aes(x = Tumor, y = Normal)) +
  geom_point(size = 1, aes(color = significance, alpha = significance)) +  # Smaller points
  scale_color_manual(values = c("Not Significant" = "gray", 
                                "Tumor Only" = "darkgreen", 
                                "Normal Only" = "blue", 
                                "Both" = "red")) +
  scale_alpha_manual(values = c("Not Significant" = 0.7, 
                                "Tumor Only" = 1, 
                                "Normal Only" = 1, 
                                "Both" = 1)) +
  # Annotate Plot
  annotate("text", x = 9, y = -2, label = "AFR T", size = 5, vjust = 1) +
  annotate("text", x = -9, y = -2, label = "EUR T", size = 5, vjust = 1) +
  annotate("text", x = -1.5, y = 6, label = "AFR N", size = 5, hjust = 0) +
  annotate("text", x = -1.5, y = -6, label = "EUR N", size = 5, hjust = 0) +
  # Add Intercepts
  # geom_vline(xintercept = 0, linetype = "solid", color = "black", size = 0.8) +
  # geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.8) +
  annotate("segment", x = -10, xend = 10, y = 0, yend = 0, color = "black", size = 0.5,
           arrow = arrow(type = "closed", ends = "both", length = unit(0.2, "cm"))) +
  annotate("segment", x = 0, xend = 0, y = -5, yend = 5, color = "black", size = 0.5,
           arrow = arrow(type = "closed", ends = "both", length = unit(0.2, "cm"))) +
  scale_x_continuous(breaks = seq(-15, 15, by = 3)) +
  scale_y_continuous(breaks = seq(-15, 15, by = 3)) +
  
  # Label points
  geom_text_repel(data = sig_both, fontface = "bold", aes(label = Gene_symbol), size = 3, max.overlaps = 30, force = 2) +
  labs(x = "Wald Statistic", y = "Wald Statistic", color = "Significance", alpha = "Significance") + 
  # Add axis ticks
  theme_minimal(base_size = 12) + 
  theme(
    axis.line = element_blank(),  
    axis.text = element_text(size = 12),
    panel.grid = element_blank(),   
    panel.border = element_blank()
  )

ggsave(plot = p, path = here(output.path), filename = "wald_AFRvsEUR_TvsN.png", width = 6, height = 4, dpi = 600, bg = "white")


## Fig 3 B Volcano Plot TCGA ####

results_df <- tcga %>% mutate(log10_pvalue = -log10(padj))
remove <- c("MTND1P23") #High cook value, driven by outlier
results_df <- results_df %>% filter(!Gene_symbol %in% remove)
label_subset <- results_df %>% filter(Gene_symbol %in% shared)

g <- ggplot(results_df, aes(x = log2FoldChange, y = log10_pvalue)) + 
  geom_point(data = subset(results_df, log2FoldChange >= 0.5 & padj < 0.05), color = "red", alpha = 0.5, aes(text = Gene_symbol)) + 
  geom_point(data = subset(results_df, log2FoldChange <= -0.5 & padj < 0.05), color = "lightblue", alpha = 0.5, aes(text = Gene_symbol)) +
  geom_hline(yintercept = -log10(.05), linetype = 3) + geom_vline(xintercept = c(-0.5,0.5), linetype = 3) + 
  geom_point(data = subset(results_df, abs(padj) >= 0.05 | abs(log2FoldChange) < 0.5 & padj <0.05), alpha = 0.05, color = "lightgray") + theme_classic() + 
  geom_text_repel(data = label_subset,fontface = "bold", max.overlaps = 30, force = 0.7, aes(x = log2FoldChange, y = log10_pvalue,  label = Gene_symbol),
                  size = 3, min.segment.length	= 0.05
  ) + 
  #coord_cartesian(xlim = c(-5, 5)) + 
  ylab(expression("-Log"[10]~Adjusted~italic(P))) + xlab(expression(Log[2])) 

ggsave(path = here(output.path), filename = "TCGA_DEA_volcano.png", plot = g, device = "png", width = 6, height = 4, units = "in", dpi = 600, bg = "white")



## GSEA Analysis for SUP FIG 3 and FIG 3 C ####

## JHU GSEA output ##
source(here('code/functions/gsea_msigdb.R'))

files <- list.files(here('output/fig3supfig3v2/Ancestry_Bin_wNormal/'), pattern = ".txt")

for(i in 1:length(files)){
  df = read_tsv(here('output/fig3supfig3v2/Ancestry_Bin_wNormal/', files[i]))
  comparison = str_extract(files[i],  "(?<=Results_).*(?=\\.txt)")
  gsea_msigdb(ranked_genes = df, gene.var = "Gene_symbol", stat = "stat", comparison = comparison)
  
}

## TCGA GSEA output ##
files <- list.files(here('output/fig3supfig3v2/TCGA_Ancestry_Tumor/'), pattern = ".txt")

for(i in 1:length(files)){
  df = read_tsv(here('output/fig3supfig3v2/TCGA_Ancestry_Tumor/', files[i]))
  comparison = str_extract(files[i],  "(?<=Results_).*(?=\\.txt)")
  gsea_msigdb(ranked_genes = df, gene.var = "Gene_symbol", stat = "stat", comparison = "GSEA_TCGA_Ancesstry_Tumor")
  
}

## Figure 3C Compare TCGA and JHU AFR_T vs EUR T GSEA ####
source(here('code/functions/gsea_multi_NES_heatmap.R'))

JHU.Hall <- read_csv(here('output/fig3supfig3v2/AFR_T_vs_EUR_T/AFR_T_vs_EUR_T_HALLMARK.csv'))
TCGA.Hall <- read_csv(here('output/fig3supfig3v2/GSEA_TCGA_Ancesstry_Tumor/GSEA_TCGA_Ancesstry_Tumor_HALLMARK.csv'))

top_genes <- list(JHU.Hall, TCGA.Hall)

JHU.Hall %>% filter(padj < 0.05)
sig <- TCGA.Hall %>% filter(padj < 0.05, NES < 0)

JHU.Hall %>% filter(padj < 0.05, pathway %in% sig$pathway)

df <- do.call(rbind, top_genes)
h <- gsea_multi_NES_heatmap(df, up.only = F)

pdf(file = here(output.path, "JHUvsTCGA_Heatmap.pdf"),
    width = 4, height = 5)
ComplexHeatmap::draw(h)
dev.off()


## Sup Fig 3 B - Comparison of GSEA of EUR T vs N and AFR T vs N in JHU####
AFR_TvsN <- read_csv(here('output/fig3supfig3v2/AFR_T_vs_AFR_N/AFR_T_vs_AFR_N_HALLMARK.csv'))
EUR_TvsN <- read_csv(here('output/fig3supfig3v2/EUR_T_vs_EUR_N/EUR_T_vs_EUR_N_HALLMARK.csv'))

top_genes <- list(AFR_TvsN, EUR_TvsN)


df <- do.call(rbind, top_genes)
h <- gsea_multi_NES_heatmap(df, up.only = F)

pdf(file = here(output.path, "TvsN_AFR_EUR_Hallmark_Heatmap.pdf"),
    width = 4, height = 5)
ComplexHeatmap::draw(h)
dev.off()

