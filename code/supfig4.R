library(matrixStats)
library(tidyverse)
library(SummarizedExperiment)
library(here)

here::i_am("code/supfig4.R")
output.path <- here("output/supfig4")
if (!dir.exists(here(output.path))) dir.create(here(output.path))

## Loading Pooled Cohort
se <- readRDS(here("output/dataprep_pooled/bc_log2tpm_pooled_se.rds"))

tpm <- assays(se)[["bc.log2.tpm"]]
keep <- rowSums(tpm) != 0


# Load IM151 Data ####


#' ### Im151 Training Dataset
#' Access to this dataset was provided though a data sharing agreement with Genentech and downloaded from the EGA (#EGAC00001001813)
#'  The data was originally described by [Motzer et al.](https://www.cell.com/cancer-cell/fulltext/S1535-6108(20)30542-0)

## Load and format
tpm.IMm151.df <- read_tsv(here("data/imm151/IMmotion151.expression.data.TPM.anon.20201106.tsv"))
tpm.IMm151 <- as.matrix(tpm.IMm151.df %>% dplyr::select(4:ncol(tpm.IMm151.df)))
rownames(tpm.IMm151) <- tpm.IMm151.df$X1

## Load Clinical Data
meta.IMm151<- read_csv(here("data/imm151/imM151_meta_clean.csv"))

## Subsetting for top 10% MAD (as described in Motzer et al.,)
tpm_rowMADs <- rowMads(tpm.IMm151)
tpm.IMm151 <- tpm.IMm151[names(sort(tpm_rowMADs, decreasing = TRUE)[1:floor(nrow(tpm.IMm151) / 10)]),] #Top 10% by MAD

## Converting GeneIDs to Gene Symbols
ImM151_genes <- select(tpm.IMm151.df, X1, symbol) %>% filter(X1 %in% rownames(tpm.IMm151))
ImM151_genes <- ImM151_genes %>% filter(!is.na(symbol))

## Identify Overlapping genes:
common <- BiocGenerics::intersect(rownames(se), ImM151_genes$symbol)
gene_id.common <- ImM151_genes %>% filter(symbol %in% common)

## Subset for overlapping genes
IM151 <- tpm.IMm151[gene_id.common$X1,]
rownames(IM151) <- gene_id.common$symbol

## Row Normalize 
IM151 <- t(scale(t(IM151
)))

##Format for RF
IM151.df <- as.data.frame(t(IM151))
IM151.df <- IM151.df %>% mutate(SampleID = colnames(IM151)) %>% dplyr::select(SampleID, everything())
IM151.df <- left_join(IM151.df, dplyr::select(meta.IMm151, SampleID, cluster)) %>% dplyr::select(SampleID, cluster, everything())

## Repair Feature Names (for Random Forest package)
fixed <- str_replace(colnames(IM151.df), "-", "_")
colnames(IM151.df) <- fixed


##Filtering out snoRNA cases 
IM151.df <- IM151.df %>% dplyr::filter(cluster != "snoRNA") #Removes 28 samples

## Reformat Cluster names for use as outcome variable in Random forest
IM151.df <- IM151.df %>% mutate(cluster = case_when(
  cluster == "Angio/Stromal" ~ "angio_stromal",
  cluster == "Angiogenic" ~ "angiogenic",
  cluster == "Complement./Ω-ox." ~ "complement_omega_oxidation",
  cluster == "Proliferative" ~ "proliferative",
  cluster == "Stromal/Proliferative" ~ "stromal_proliferative",
  cluster == "T-eff/Proliferative" ~ "teff_proliferative",
  cluster == "snoRNA" ~ "snoRNA"
))
IM151.df$cluster <- factor(IM151.df$cluster, levels = c("angiogenic", "angio_stromal", "complement_omega_oxidation", "proliferative", "stromal_proliferative", "teff_proliferative"))

## Saving Training/validation set
write_csv(IM151.df, here(output.path, "IM151_rf_dataset.csv"))

# Format JHU cohort ####

## Subset for shared genes
se <- se[gene_id.common$symbol,]

## Log Transform and row normalize
log.tpm <- assays(se)[["bc.log2.tpm"]] 
log.tpm <- t(scale(t(log.tpm)))

## Re-format as dataframe
pooled.df <- as.data.frame(t(log.tpm))
pooled.df <- pooled.df %>% mutate(SampleID = colnames(log.tpm)) %>% dplyr::select(SampleID, everything())

## Repair feature(gene) names 
fixed <- str_replace(colnames(pooled.df), "-", "_")
colnames(pooled.df) <- fixed

#Save
write_csv(pooled.df, here(output.path, "pooled_rf_dataset.csv"))


# Train RF Model ####
library(caret)
library(ranger)
library(rsample)
library(doParallel)


set.seed(1353)
split <- initial_split(IM151.df, prop = 0.8, 
                       strata = "cluster")
IM151_train  <- training(split)
IM151_test   <- testing(split)

IM151_test %>% write_csv(here(output.path, "IM151_test_dataset.csv"))
IM151_train %>% write_csv(here(output.path, "IM151_train_dataset.csv"))

## Cross Validation training settings
cv_control <- trainControl(method = "cv", 
                           number = 5, 
                           classProbs = TRUE,
                           summaryFunction = multiClassSummary, 
                           savePredictions = "final", # Save predictions for the best model
                           sampling = "up", #Upsampling since the data is unblanced
                           verboseIter = TRUE,
                           returnResamp = "all",
                           allowParallel = TRUE)
## Tune Paramters: 

hyper_grid <- expand.grid(mtry = c(25, 50, 100, 120, 140, 250, 600),
                          min.node.size = c(1, 5, 10,25,50), 
                          splitrule = "gini")

cl <- makePSOCKcluster(4)
registerDoParallel(cl)
model_tuning <- train(cluster ~ .,
                      data = IM151_train,
                      method = "ranger",
                      trControl = cv_control,
                      tuneGrid = hyper_grid,
                      num.trees = 2700,
                      metric = "Accuracy")
stopCluster(cl)


tuning_results <- model_tuning$results
tuning_results %>% write_csv(here(output.path, "random_forest_tuning_results.csv"))



tuning_results <- read_csv(here("output/supfig4/random_forest_tuning_results.csv"))


# Visualizing Tunirng Results ####
tuning_results <- tuning_results %>% mutate(node = tuning_results$min.node.size,
                                            mtry_factor = tuning_results$mtry)


g <- ggplot(tuning_results, aes(x = min.node.size, y = Accuracy)) +
  geom_point(aes(color = factor(mtry)), alpha = 0.6) +
  geom_line(aes(group = factor(mtry), color = factor(mtry))) +
  scale_color_brewer(palette = "Set1", name = "mtry") +
  theme_minimal() +
  labs(
    title = "mtry and min.node.size",
    x = "min.node.size",
    y = "Accuracy"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.position = "bottom"
  )


ggsave(filename = "rf_tuning.pdf", path = output.path, plot = g, device = "pdf", width = 5, height = 3.5, unit = "in")


# Training RF model on optimal tuning results ####
IM151_train <- read_csv(here("output/supfig4/IM151_train_dataset.csv"))

IM151_train$cluster <- factor(IM151_train$cluster)

set.seed(12345)
fit <- ranger(
  formula         = cluster ~ ., 
  data            = IM151_train, 
  num.trees       = 2000,
  mtry            = 50,
  min.node.size   = 10,
  replace = TRUE,
  verbose = TRUE,
  num.threads = 6
)

fit %>% saveRDS(here(output.path, "rf_model.rds"))


library(ranger)
library(caret)
library(tidyverse)
library(here)

# Assessing performance on Validation Cohort ####
fit <- readRDS(here('output/supfig4/rf_model.rds'))
IM151_test <- read_csv(here("output/supfig4/IM151_test_dataset.csv"))
IM151_test$cluster <- factor(IM151_test$cluster)
pred.IM151 <- predict(fit, data = IM151_test)
predictions <- tibble(SampleID = IM151_test$SampleID, cluster = IM151_test$cluster, prediction = pred.IM151$predictions)



#Confusion Matrix
confusion_matrix <- confusionMatrix(pred.IM151$predictions, IM151_test$cluster)

#Overall Performance
confusion_matrix$overall

#Performance by molecular class
df <-  as.data.frame(confusion_matrix$byClass)
df$cluster <- rownames(df)
write_csv(df, here(output.path, "rf_class_performance.csv"))

# Applying to JHU Cohort ####
jhu.df <- read_csv(here('output/supfig4/pooled_rf_dataset.csv'))

jhu_predictions <- predict(fit, data = jhu.df)
results <- tibble(cluster = jhu_predictions$predictions, ptID = jhu.df$SampleID)


se <- readRDS(here("output/dataprep_pooled/bc_log2tpm_pooled_se.rds"))
write_csv(results, here(output.path, "IMM151_predications.csv"))

sessionInfo()
