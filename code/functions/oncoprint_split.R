#' @title Heatmap Oncoprint Default
#' @param input.matrix Matrix to be plotted. Must contain colnames equal to all columns of matrices 
#' @param input.meta Dataframe which provides top annotation information
#' @param goi Genes to subset matrix on 
#' @param split_var metadata variable to split by
#' @param ... Additional parameters
#' @return  Oncoprint


oncoprint_split <- function(input.matrix, input.meta, goi, split_var, top_anno, ...){
  source(here('code/functions/oncoprint_default.R'))
  source(here('code/functions/heatmap_top_annotation_colors.R'))
  source(here('code/functions/stats_group_frequencies_and_pvals.R'))
  
  
  #  Top Annotation
  subtype.col <- c("Angiogenic" = "#6d1f82",
                   "Angio/Stromal" =  "#fa170f",
                   "Stromal/Proliferative" = "gray3",
                   "Proliferative" = "#faa80f", 
                   "T-eff/Proliferative" = "#1f8248", 
                   "Complement/Omega-ox." = "#539bed")
  
  
  source(here('code/functions/heatmap_top_annotation_colors.R'))
  
  
  input.matrix <- input.matrix[goi,]
  stats <- stats_group_frequencies_and_pvals(input.matrix = input.matrix, meta_data = input.meta, comparison_variable = split_var)
  
  split <- sort(unique(input.meta[[split_var]]))
  
  
  annotation_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "bold"), 
                                 grid_height = unit(2, "mm"),
                                 grid_width = unit(2, "mm"), 
                                 direction = "horizontal",
                                 labels_gp = gpar(fontsize = 8), 
                                 annotation_name_align = TRUE)
  oncoprint_list <- list()
  
  for(i in split){
    print(i)
    split_data <- input.meta %>% filter(!!sym(split_var) == i)
    split_mat <- input.matrix[,split_data$sampleID]
    split_anno <- split_data %>% select(-sampleID)
    
    if(i == split[1]){
      freq_var <- paste0(i, "_freq")
      formatted_percentages <- paste0(round(stats[[freq_var]] * 100), "%")
      r_anno = rowAnnotation(annotation_name_gp = gpar(fontsize = 8),
                             show_annotation_name = FALSE,
                             "Pct" = anno_text(
                               formatted_percentages, 
                               gp = gpar(fontsize = 8, fontface = "bold"),
                             ))
      l_anno = rowAnnotation(annotation_name_gp = gpar(fontsize = 8),
                             show_annotation_name = FALSE,
                             "Gene" = anno_text(
                               stats$gene, 
                               gp = gpar(fontsize =8, fontface = "bold"),
                             ))
      
      top_anno <- HeatmapAnnotation(
        height = unit(2, "cm"),
        simple_anno_size = unit(0.2, "cm"),
        annotation_legend_param = annotation_legend_param,
        show_legend = TRUE,
        show_annotation_name = FALSE,
        annotation_name_gp = gpar(fontsize = 8),
        AFR = split_anno$AFR,
        EUR = split_anno$EUR,
        Subtype = split_anno$cluster_label,
        col = list(AFR = AFR.col, EUR = EUR.col,
                   Subtype = subtype.col))
    } else if(i != split[length(split)]){
      freq_var <- paste0(i, "_freq")
      formatted_percentages <- paste0(round(stats[[freq_var]] * 100), "%")
      r_anno = rowAnnotation(annotation_name_gp = gpar(fontsize =8),
                             show_annotation_name = FALSE,
                             "Pct" = anno_text(
                               formatted_percentages, 
                               gp = gpar(fontsize = 8, fontface = "bold"),
                             ))
      l_anno = NULL
      
      top_anno <- HeatmapAnnotation(
        height = unit(2, "cm"),
        simple_anno_size = unit(0.2, "cm"),
        annotation_legend_param = annotation_legend_param,
        show_legend = TRUE,
        show_annotation_name = FALSE,
        annotation_name_gp = gpar(fontsize = 8),
        AFR = split_anno$AFR,
        EUR = split_anno$EUR,
        Subtype = split_anno$cluster_label,
        col = list(AFR = AFR.col, EUR = EUR.col,
                   Subtype = subtype.col))
    } else{
      freq_var <- paste0(i, "_freq")
      formatted_percentages <- paste0(round(stats[[freq_var]] * 100), "%")
      r_anno = rowAnnotation(annotation_name_gp = gpar(fontsize = 8),
                             show_annotation_name = FALSE,
                             "Pct" = anno_text(
                               formatted_percentages, 
                               gp = gpar(fontsize = 8, fontface = "bold")),
                               "Pval" = anno_text(
                                 stats$pval_lab, 
                                 gp = gpar(fontsize = 8, fontface = "bold"))
                               )
      l_anno = NULL
      top_anno = HeatmapAnnotation(
        height = unit(2, "cm"),
        simple_anno_size = unit(0.2, "cm"),
        annotation_legend_param = annotation_legend_param,
        show_legend = TRUE,
        show_annotation_name = FALSE,
        annotation_name_gp = gpar(fontsize = 8),
        AFR = split_anno$AFR,
        EUR = split_anno$EUR,
        Subtype = split_anno$cluster_label,
        col = list(AFR = AFR.col, EUR = EUR.col,
                   Subtype = subtype.col))
    }
    oncoprint_list[[i]] <- 
      oncoprint_default(input.matrix = split_mat, 
                                input.meta = split_data, 
                                show_pct = F,
                                show_row_names = F,
                                top_anno = top_anno,
                                left_annotation = l_anno,
                                right_annotation = r_anno,
                                row_names_side = "left")
  }
  
  combined <- Reduce("+", oncoprint_list)
  return(combined)
}

# oncoprint_default(input.matrix = split_mat, 
#                   input.meta = split_data, 
#                   show_pct = F,
#                   show_row_names = F,
#                   top_anno = top_anno,
#                   left_annotation = l_anno,
#                   right_annotation = r_anno,
#                   row_names_side = "left")
