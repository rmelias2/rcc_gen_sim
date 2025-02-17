#'Oncoprint_default

#' @param input.matrix Matrix to be plotted. Must contain colnames equal to all columns of matrices 
#' @param input.meta Dataframe which provides top annotation information
#' @param title Title for heatmap legend
#' @param anno.legend Annotation legend paramaters. Default Filled below
#' @param hm.legend description
#' @param top_anno description

oncoprint_default <- function(input.matrix, input.meta, title, anno.legend = "default", 
                              hm.legend = "default", top_anno = "default", ...){
  
  required_packages <- c("ComplexHeatmap", "RColorBrewer", "circlize")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("Package", pkg, "is required but not installed. Please install it using install.packages()."))
    } else {
      library(pkg, character.only = TRUE)
    }
  }
  
  source(here('code/functions/heatmap_top_annotation_colors.R'))
  source(here('code/functions/oncoprint_variant_annotations.R'))
  
  #Assign Variant Labels
  unique_alterations <- function(matrix) {
    flattened_vector <- na.omit(c(matrix))
    split_values <- strsplit(flattened_vector, split = ",", fixed = TRUE)
    unique_values <- unique(unlist(split_values))
    return(unique_values)
  }
  
  variant_labels <- names(graphics.final)[2:length(names(graphics.final))]
  variant_labels <- variant_labels[variant_labels %in% unique_alterations(input.matrix)]
  
  # Verify Input meta Data matches Input Matrix
  if(!is.null(input.meta)){
    coldata = input.meta
    
    if(all(colnames(input.matrix) != coldata$sampleID)){
      stop("Input Matrix and Input Metadata files do not match")
    }
  }
  
  #Set Annotation Legend
  if(anno.legend == "default"){
    annotation_legend_param = list(title_gp = gpar(fontsize = 4, fontface = "bold"), 
                                   grid_height = unit(2, "mm"),
                                   grid_width = unit(2, "mm"), 
                                   ncol = 2,
                                   direction = "horizontal",
                                   labels_gp = gpar(fontsize = 4), 
                                   annotation_name_align = TRUE)
  }
  
  if(hm.legend == "default"){
    hm.legend = list(title = "Variant",
                     title_gp = gpar(fontsize = 4, fontface = "bold"), 
                     grid_height = unit(2, "mm"),
                     grid_width = unit(2, "mm"),  
                     labels_gp = gpar(fontsize = 4),
                     labels = variant_labels
)
  }
    
  if (is.character(top_anno) && top_anno == "default") {
    top_anno <- HeatmapAnnotation(
      # height = unit(2, "cm"), # Uncomment if needed
      simple_anno_size = unit(0.2, "cm"),
      annotation_legend_param = annotation_legend_param,
      show_legend = TRUE,
      annotation_name_gp = gpar(fontsize = 6),
      AFR = coldata$AFR,
      EUR = coldata$EUR,
      col = list(
        AFR = AFR.col, 
        EUR = EUR.col)
      )
  }
  #Set Arguments for oncoprint
  hm_args <- list(
    mat = input.matrix,
    alter_fun = graphics.final,
    heatmap_legend_param = hm.legend,
    column_names_gp = gpar(fontsize = 4, fontface = "bold"),
    row_names_gp = gpar(fontsize = 6, fontface = "bold"),
    ...
  )
    
    if(!is.null(input.meta)){
      hm_args$top_annotation = top_anno
    }
    
    p <-  do.call("oncoPrint", hm_args)
    return(p)
   
}
    
    
    
    
    