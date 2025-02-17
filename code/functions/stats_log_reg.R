#' stats_log_reg
#'
#' @param input.matrix Gene/Feature x Sample. Un-altered = NA, altered = !NA
#' @param meta_data LONG data format
#' @param comparison_variable Variable containing the comparison value names (in long format)
#' @param value_variable Variable containing the corresponding values
#' @param sampleID_name 
#'
#' @return
#' @export
#'
#' @examples
stats_log_reg <- function(input.matrix, meta_data, comparison_variable, sampleID_name = "sampleID", value_variable = "Value") {
  # Step 1: Validate Inputs
  if (!is.matrix(input.matrix)) stop("input.matrix should be a matrix.")
  if (!is.data.frame(meta_data)) stop("meta_data should be a data frame.")
  if (!comparison_variable %in% names(meta_data)) stop("comparison_variable not found in meta_data.")
  
  # Step 2: Prepare Output Template
  unique_patterns <- unique(meta_data[[comparison_variable]])
  output_columns <- c("Gene", 
                      sapply(unique_patterns, function(x) c(paste0(x, "_OR"), paste0(x, "_p"), paste0(x, "_adj.p"))))
  output_columns <- unlist(output_columns)
  
  
  output.df <- data.frame(matrix(ncol = length(output_columns), nrow = nrow(input.matrix)))
  colnames(output.df) <- output_columns
  output.df$Gene <- rownames(input.matrix)

  # Step 3: Loop Over Each Gene and Pattern
  for (i in seq_len(nrow(input.matrix))) {
    gene_mutation_status <- input.matrix[i,]
    
    for (pattern in unique_patterns) {
      # Subset metadata for the current pattern
      meta_subset <- meta_data[meta_data[[comparison_variable]] == pattern, ]
      
      # Merge mutation data with the meta data
      # NA values in mutation data become "FALSE [WT], and !NA beome TRUE [Mut]
      data_for_glm <- data.frame(sampleID = colnames(input.matrix),
                                 Mutation_Status = !is.na(gene_mutation_status))
      
      data_for_glm <- merge(data_for_glm, meta_subset, by.x = "sampleID", by.y = sampleID_name)
      
      # Fit binomial GLM
      # This Requires "Value"
      fit <- try(glm(Mutation_Status ~ get(value_variable), data = data_for_glm, family = binomial()), silent = TRUE)
      
      if (inherits(fit, "try-error")) {
        warning(paste("GLM failed for gene:", rownames(input.matrix)[i], "and pattern:", pattern))
        output.df[i, paste0(pattern, "_OR")] <- NA
        output.df[i, paste0(pattern, "_p")] <- NA
      } else {
        # Extract OR, p-value
        coef_summary <- summary(fit)$coefficients
        or_value <- exp(coef_summary[2, 1])
        p_value <- coef_summary[2, 4]
        
        output.df[i, paste0(pattern, "_OR")] <- or_value
        output.df[i, paste0(pattern, "_p")] <- p_value
      }
    }
    # Step 4: Adjust p-values for multiple testing (per gene)
    p_values <- unlist(output.df[i, grep("_p$", colnames(output.df))])
    adjusted_p_values <- p.adjust(p_values, method = "BH")
    output.df[i, grep("_adj.p$", colnames(output.df))] <- adjusted_p_values
  }
  return(output.df)
}
