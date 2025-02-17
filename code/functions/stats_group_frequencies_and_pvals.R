#' Function to Calculate Frequency of binary (Gain or Not, Mut or Not) accross groups 
#' @param input.matrix matrix containing NA, Missing and 0 values (Negative) or anything else (alteration)
#' @param meta_data Metadata file which contains the samples names and the grouping categories
#' @param comparison_variable Column name with the comparison variable. Can have as many levels is wanted. 
#' @param sampleID_name Default "sampleID". Column name which contains the sample IDs which match the column names of the oncoprint matrix

stats_group_frequencies_and_pvals <- function(input.matrix, meta_data, comparison_variable, sampleID_name = "sampleID"){
  if(length(colnames(input.matrix)) != length(meta_data[[sampleID_name]])) {
    stop("Input Matrix and Input Metadata files do not match in length.")
  } else if(!all(colnames(input.matrix) == meta_data[[sampleID_name]])) {
    # If lengths are the same, check if all values are equal
    stop("Input Matrix and Input Metadata files do not match in values.")
  }
  
  ## Prepare input
  comparison_variable <- sym(comparison_variable)
  sampleID <- sym(sampleID_name)
  df <- as.data.frame(t(input.matrix))
  genes  <- rownames(input.matrix)
  colnames(df) <- genes
  df[df == "" | is.na(df) | df == 0 | df == FALSE] <- "Wt"  
  df[df!= "Wt"] <- "Mut" #Simplifies the above matrix which has the type of mutation included
  df <- df %>%
    mutate(!!sampleID := rownames(df))
  meta_data[[sampleID]] <- as.character(meta_data[[sampleID]])
  df[[sampleID]] <- as.character(df[[sampleID]])
  df <- left_join(df,select(meta_data, !!sampleID,!!comparison_variable)) %>% select(!!sampleID, !!comparison_variable, everything())
  
  ## Prepare dataframe to collect statistics 
  unique_groups <- unique(df[[comparison_variable]])
  if(length(unique_groups) == 2) {
    print("The comparison variable has two groups, performing the fisher exact test")
    group1 <- unique_groups[1]
    group2 <- unique_groups[2]
    group_specific_colnames <- c("count", "freq")
    colnames_result <- c(
      paste(group1, group_specific_colnames, sep = "_"),
      paste(group2, group_specific_colnames, sep = "_")
    )
    colnames_result <- c(colnames_result, "OR", "p_val","CI_low", "CI_high")
    results <- tibble(gene = genes)
    results[colnames_result] <- 0 
    
    for(i in seq_along(genes)){
      cont_tbl <- table(df[, genes[i]], df[[comparison_variable]])
      if (nrow(cont_tbl) == 1) {
        stop(paste0("No alterations detected in feature: ", genes[i]))
      }
      cont_df <- as.data.frame(as.table(cont_tbl))  # Convert cont_tbl to a data frame
      present_df <- cont_df %>% filter(Var1 == "Mut")
      group1_count <- sum(present_df$Freq[present_df$Var2 == group1])
      group1_freq <- group1_count / sum(cont_df$Freq[cont_df$Var2 == group1])
      
      group2_count <- sum(present_df$Freq[present_df$Var2 ==  group2])
      group2_freq <- group2_count /  sum(cont_df$Freq[cont_df$Var2 == group2])
      
      results[[paste(group1, "count", sep = "_")]][i] <- group1_count
      results[[paste(group1, "freq", sep = "_")]][i] <- group1_freq
      results[[paste(group2, "count", sep = "_")]][i] <- group2_count
      results[[paste(group2, "freq", sep = "_")]][i] <- group2_freq
      
      res <- fisher.test(cont_tbl) 
      results$OR[i] <- res$estimate
      results$CI_low[i] <- res$conf.int[1]
      results$CI_high[i] <- res$conf.int[2]
      results$p_val[i] <- res$p.value
    }
    results$adj.p <- p.adjust(results$p_val, method = "fdr")
    results[[paste(group1, "freq", sep = "_")]] <- round(results[[paste(group1, "freq", sep = "_")]], 3)
    results[[paste(group2, "freq", sep = "_")]] <- round(results[[paste(group2, "freq", sep = "_")]], 3)
    results <- results %>% mutate(p_val = round(p_val, digits = 3),
                                  pval_lab = if_else(p_val < 0.0001, "     <0.0001", paste0("    ",p_val)))
  }  
  else if(length(unique_groups) > 2) {
    print(paste0("There are ", length(unique_groups), " groups. Returning only Pvalue"))
    group_specific_colnames <- c("count", "freq")
    colnames_result <- unlist(lapply(unique_groups, function(g) paste(g, group_specific_colnames, sep = "_")))
    results <- tibble(gene = genes)
    results[colnames_result] <- 0  # Initialize columns for counts and frequencies
    
    for (i in seq_along(genes)) {
      cont_tbl <- table(df[, genes[i]], df[[as.character(comparison_variable)]])
      if (nrow(cont_tbl) == 1) {
        stop(paste0("No alterations detected in feature: ", genes[i]))
      }
      
      for (group in unique_groups) {
        group_mut_count <- sum(cont_tbl["Mut", group, drop = FALSE], na.rm = TRUE)
        group_total_count <- sum(cont_tbl[, group, drop = FALSE], na.rm = TRUE)
        group_freq <- group_mut_count / group_total_count
        
        results[[paste(group, "count", sep = "_")]][i] <- group_mut_count
        results[[paste(group, "freq", sep = "_")]][i] <- group_freq
      }
      
      # Results returned depend on the number of unique groups
      if (length(unique_groups) == 2) {
        res <- fisher.test(cont_tbl)
        results$OR[i] <- res$estimate
        results$CI_low[i] <- res$conf.int[1]
        results$CI_high[i] <- res$conf.int[2]
        results$p_val[i] <- res$p.value
      } else {
        res <- fisher.test(cont_tbl, simulate.p.value=TRUE)
        results$p_val[i] <- res$p.value
      }
    }
    results$adj.p <- p.adjust(results$p_val, method = "fdr")
    results <- results %>%
      mutate(across(ends_with("freq"), round, 3),
             p_val = round(p_val, 3),
             pval_lab = if_else(p_val < 0.0001, "     <0.0001", sprintf("    %.4f", p_val)))
    
    return(results)
  }
}

