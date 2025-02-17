library(here)
library(tidyverse)

here::i_am("code/fig2.R")
output.path <- here("output/fig2")
if (!dir.exists(here(output.path))) dir.create(here(output.path))


#Ancestry Distribution Fig2A ####
meta <- read_csv(here('output/table2/pooled_cohort.csv'))
plot <- meta %>% filter(!is.na(Race))

# Visualize Ancestry Distribution: 

g <- ggplot(plot, aes(x = AFR, fill = Race)) +
  scale_fill_manual(values = c("B" = "red", "W" = "blue")) + 
  geom_density(alpha = 0.5, position = "identity", color = "black", binwidth = 0.05) + 
  labs(
    x = "% AFR Ancestry",
    y = "Count") +
  theme_minimal() +
  theme(legend.position = "top", 
        legend.title = element_blank(),
        axis.title = element_text(size = 8), 
        axis.text = element_text(size = 8), 
        legend.text = element_text(size = 8)) 

ggsave(path = output.path, filename = "AFR_Density_Plot.pdf", plot = g, device = "pdf", width = 3, height = 2.5, units = "in")



# Group Level Comparisons ####

# Identify somatic mutation frequency differences by group  (2B) ####
# Load Functions
source(here('code/functions/format_maf_to_matrix.R'))
source(here('code/functions/stats_group_frequencies_and_pvals.R'))
source(here('code/functions/ggplot_alteration_frequency_scatter.R'))
source(here('code/functions/stats_log_reg.R'))


## Function for preparing statistical comparison ####
stats_group_data_prep <- function(input.maf = NULL, input.matrix = NULL, input.matrix.cuttoff = NULL, meta_data, meta_sampleID = "sampleID", maf_sampleID = "Tumor_Sample_Barcode",
                                  min_alterations = 3, comparison_variable){
  samples <- unique(meta_data[[meta_sampleID]])
  
  #Extract input.matrix.cuttoff
  if (!is.null(input.matrix.cuttoff)) {
    match <- regexpr("<=|>=|<|>|==|!=", input.matrix.cuttoff)
    operator <- regmatches(input.matrix.cuttoff, match)
    cutoff <- as.numeric(sub(operator, "", input.matrix.cuttoff))
  }
  if(!is.null(input.maf)){
    maf.subset <- maf %>% filter(maf[[maf_sampleID]] %in% meta_data[[meta_sampleID]])
    a <- maf.subset %>% group_by(Hugo_Symbol, maf[[maf_sampleID]]) %>% slice_head() %>% 
      ungroup() %>%  dplyr::count(Hugo_Symbol) %>% arrange(desc(n))
    a <- a %>% filter(n >= min_alterations)
    goi <- a$Hugo_Symbol
    
    
    gene_level_mat <- format_maf_to_matrix(maf.df = maf.subset, genes = goi, sampleID = maf_sampleID)
    full_mat <- annotate_no_call(gene_level_mat, samples)
    
  } else if(!is.null(input.matrix)){
    # Apply the cutoff condition based on the extracted operator
    cn.mat <- input.matrix
    if (operator == ">") {
      cn.mat[cn.mat > cutoff] <- "del"
    } else if (operator == "<") {
      cn.mat[cn.mat < cutoff] <- "del"
    } else if (operator == ">=") {
      cn.mat[cn.mat >= cutoff] <- "del"
    } else if (operator == "<=") {
      cn.mat[cn.mat <= cutoff] <- "del"
    } else {
      stop("Unsupported operator. Use one of '>', '<', '>=', '<='.")
    }
    
    cn.mat[cn.mat != "del"] <- NA
    full_mat <- cn.mat
    
  }
  
  #Subset mat for samples of interest
  full_mat <- full_mat[,samples]
  
  #Exclude rows with less than min alterations
  goi <- rowSums(!is.na(full_mat)) %>% sort(decreasing = TRUE) 
  goi <- goi[goi > min_alterations]
  comparison.mat <- full_mat[names(goi),]
  
  stats <- stats_group_frequencies_and_pvals(input.matrix = comparison.mat, meta_data = meta_data, comparison_variable = comparison_variable, sampleID_name = meta_sampleID)
  return(stats)
}


## SOMATIC MUTAITON COMPARISON ####
### Group Level
maf <- read_tsv(here('output/dataprep_pooled/pooled_maf.tsv'))

comparison.df <- meta %>% filter(ptID %in% maf$ptID)
stats <- stats_group_data_prep(input.maf = maf, 
                               input.matrix = NULL,
                               meta_data = comparison.df, 
                               comparison_variable = "ancestry_group", 
                               meta_sampleID = "ptID",min_alterations = 3, maf_sampleID = "ptID")
write_csv(stats, here(output.path, "stats_muts_AFRvsEUR.csv"))

labels <- c("VHL", "PBRM1", "BAP1", "SETD2", "TTN")

stats <- stats %>% mutate(FDR = if_else(adj.p < 0.1, "<0.1", "ns"))
g <- ggplot_alteration_frequency_scatter(plot = stats, x_val = "EUR_freq", y_val = "AFR_freq", 
                                         color_var = "FDR",labels = labels, color = "FDR",
                                         x = "EUR", y = "AFR", xlim = 0.6, ylim = 0.6)

ggsave(plot = g, filename = "mutations_EURvsAFR_pooled.pdf", path = here(output.path), width = 4, height = 3, device = "pdf")


## Cytoband level Comparisons (2C)####
### Group level comparison
cytoband <- read_rds(here('output/dataprep_pooled/pooled_cns_se.rds'))
cytoband.mat <- assays(cytoband)[[1]]
comparison.df <- meta %>% filter(ptID %in% colnames(cytoband.mat))

gains <- cytoband.mat 
gains[gains >= 0.3] <- "alt"
del <- cytoband.mat 
del[del <= -0.3] <- "alt"


rownames(gains) <- paste0("gain ",rownames(gains))
rownames(del) <- paste0("del ",rownames(del))

cn.mat <- rbind(gains, del)
cn.mat[cn.mat != "alt"] <- NA
goi <- rowSums(!is.na(cn.mat)) %>% sort(decreasing = TRUE) 
goi <- goi[goi > 3]
cn.mat <- cn.mat[names(goi),]


stats <- stats_group_frequencies_and_pvals(input.matrix = cn.mat, meta_data = comparison.df, comparison_variable = "ancestry_group", sampleID_name = "ptID")

write_csv(stats, here(output.path, "stats_CN_AFRvsEUR.csv"))

labels <- c("del chr14q21.3", "del chr3p21.1", "del chr9p21.3", "del chr3p21.33", "gain chr5q35.1")

stats <- stats %>% mutate(FDR = if_else(adj.p < 0.1, "<0.1", "ns"))
g <- ggplot_alteration_frequency_scatter(plot = stats, x_val = "EUR_freq", y_val = "AFR_freq", 
                                         color_var = "FDR", labels = labels, color = "FDR",
                                         x = "EUR", y = "AFR", xlim = 1, ylim = 1)

ggsave(plot = g, filename = "CN_EURvsAFR_pooled.pdf", path = here(output.path), width = 4, height = 3, device = "pdf")

### Figure 2D New (VHL)
afr.data <- comparison.df %>% filter(Race == "B")
afr.data <- afr.data %>% select(ptID, AFR)

vhl.status = tibble(ptID = colnames(gene_level_mat), vhl = if_else(is.na(gene_level_mat["VHL",]), 0, 1), pbrm1 = if_else(is.na(gene_level_mat["PBRM1",]), 0, 1))
chr3p.status = tibble(ptID = colnames(gene_level_mat), chr3p = if_else(is.na(cn.mat["del chr3p21.1",]), 0, 1))


afr.data <- left_join(afr.data, vhl.status)
afr.data <- left_join(afr.data, chr3p.status)

g <- ggplot(afr.data, aes(x = AFR, y = vhl)) +
  geom_point()  

model <- glm(vhl ~ AFR, data = afr.data, family = binomial)
summary(model)
sessionInfo()

ggplot(afr.data, aes(x = AFR, y = chr3p)) +
  geom_point(position = position_jitter(height = 0.02)) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"))


#Figure 2D
afr.data <- comparison.df %>% filter(ancestry_group == "AFR")
afr.data <- afr.data %>% select(ptID, AFR)
afr.data$group = "AFR"

cn.mat <- rbind(gains, del)
cn.mat <- cn.mat[,afr.data$ptID]
cn.mat[cn.mat != "alt"] <- NA
goi <- rowSums(!is.na(cn.mat)) %>% sort(decreasing = TRUE) 
goi <- goi[goi > 3]
cn.mat <- cn.mat[names(goi),]

maf.subset <- maf %>% filter(ptID %in% afr.data$ptID)
a <- maf.subset %>% group_by(Hugo_Symbol, ptID) %>% slice_head() %>% 
  ungroup() %>%  dplyr::count(Hugo_Symbol) %>% arrange(desc(n))
a <- a %>% filter(n >= 3)
goi <- a$Hugo_Symbol

gene_level_mat <- format_maf_to_matrix(maf.df = maf.subset, genes = goi, sampleID = "ptID")
full_mat <- annotate_no_call(x = gene_level_mat, samples = afr.data$ptID)

vhl.status = tibble(ptID = colnames(full_mat), vhl = if_else(is.na(full_mat["VHL",]), 0, 1))

afr.data <- left_join(afr.data, vhl.status)

model <- glm(vhl ~ AFR, data = afr.data, family = binomial)
summary(model)
p_val <- summary(model)$coefficients["AFR","Pr(>|z|)"]

p_label <- if (p_val < 0.05) {
  paste0("p = ", format(p_val, digits = 2))  # round or format as you wish
} else {
  "n.s."
}

g <- ggplot(afr.data, aes(x = AFR, y = vhl)) +
  geom_point(
    position = position_jitter(height = 0.02), 
    alpha = 0.6, 
    size = 2, 
    color = "steelblue"
  ) +
  geom_smooth(
    method = "glm",
    method.args = list(family = "binomial"),
    color = "firebrick",
    se = FALSE
  ) +
  scale_y_continuous(
    breaks = c(0, 1),
    labels = c("WT", "Mut")
  ) +
  labs(
    x = "AFR %",
    y = NULL
  ) +

  annotate("text", x = 0.8, y = 0.9, 
           label = p_label, 
           color = "black", 
           size = 5) +
  theme_minimal() 
ggsave(plot = g, filename = "VHLvsAFR.pdf", path = here(output.path), device = "pdf", width = 3, height = 3)
