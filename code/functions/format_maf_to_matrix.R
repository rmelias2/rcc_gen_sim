#' Function to generate matrix of fusion and over expression calls
#' @param maf.df data frame with mutation information. Must contain variable "Variant_Classification" 
#' which contains the possible classes of coding domain sequences that are specified in the classifier function
#' @param genes vector of genes for fusion matrix 
#' @param sampleID Default `Tumor_Sample_Barcode`. Column name with sample IDs
#' @param classifier Default "Simple". Options are "Custom" and "Detailed". If Custom, must supply a list of new classifiers to the next argument
#' @param custom.classifier Default NULL. Must provide a list of new classifers
#' @param gene_lvl_cn Default NULL. Otherwise dataframe of gene level copy number changes. Must contain column "gene", "log2"
#' @param gene_lvl_cn_sampleID default "sampleID" column which names the gene level sampleID

format_maf_to_matrix <- function(maf.df, sampleID = "Tumor_Sample_Barcode", geneID = "Hugo_Symbol", genes = NULL, 
                              classifier = c("Custom", "Simple", "Detailed")[2], custom.classifier = NULL,
                              gene_lvl_cn = NULL, gene_lvl_cn_sampleID = "sampleID", gene_lvl_cn_geneID = "gene", 
                              deep_del_cutoff = -1, hetero_del_cuttoff = -0.2) {
    
  # Recode mutation classifications 
  if(classifier == "Simple"){
      classification_list = list(
        "Frameshift" = c("Frame_Shift_Del", "Frame_Shift_Ins"),
        "Indel" = c("In_Frame_Del", "In_Frame_Ins"),
        "SNV" = c("Missense_Mutation", "Nonsense_Mutation","Nonstop_Mutation", "Translation Start Site"),
        "Splice Region" = c("Splice_Site", "Splice_Region"))
    } 
  if(classifier == "Detailed"){
    classification_list = list(
      "Frameshift" = c("Frame_Shift_Del", "Frame_Shift_Ins"),
      "In-frame Indel" = c("In_Frame_Del", "In_Frame_Ins"),
      "Missense" = "Missense_Mutation",
      "Nonsense" = "Nonsense_Mutation",
      "Nonstop" = "Nonstop_Mutation",
      "Splice Region" = c("Splice_Site", "Splice_Region"),
      "Translation Start Site" = "Translation_Start_Site")
  }
  if(classifier == "Custom"){
    if(is.null(custom.classifier)){
      stop("Must define classifiers if setting classifier to custom")
    }
    classification_list = custom.classifier
  }                            
    # Initialize a column for variant labels
  maf.df$variant <- NA
  if(!is.null(genes)){
    #Subset to specific genes
  maf.df <- maf.df %>% filter(.data[[geneID]] %in% genes)
  }
    
    # Iterate over the classification list and update the variant_labels column
    for (classification in names(classification_list)) {
      labels <- classification_list[[classification]]
      for (label in labels) {
        maf.df$variant[maf.df$Variant_Classification %in% label] <- classification
      }
    }
    mut_mat <-  maf.df %>% group_by(!!sym(sampleID), !!sym(geneID)) %>% 
      summarise(variant = case_when(length(variant) > 1 ~ "Multihit", 
                                           length(variant) == 1 ~ variant))
    mut_mat <- unique(mut_mat)
    
    
    # Add gene level copy number information 
    if(!is.null(gene_lvl_cn)){
      gene_lvl_deep_del <- gene_lvl_cn %>% filter(log2 <= deep_del_cutoff, gene %in%  genes)
      gene_lvl_deep_del <- gene_lvl_deep_del %>% select(gene, sampleID)
      gene_lvl_deep_del$variant <- "Deep Deletion"
      gene_level_deletion_df <- gene_lvl_deep_del
      
      if(!is.null(hetero_del_cuttoff)){
      gene_lvl_shallow_del <- gene_lvl_cn %>% filter(log2 <= hetero_del_cuttoff & log2 > deep_del_cutoff, gene %in%  genes)
      gene_lvl_shallow_del <- gene_lvl_shallow_del %>% select(gene, sampleID)
      gene_lvl_shallow_del$variant <- "Shallow Deletion"
      gene_level_deletion_df <- bind_rows(gene_lvl_shallow_del, gene_lvl_deep_del)
      }
     
      gene_level_deletion_df <- gene_level_deletion_df %>%
        select(setNames(all_of(gene_lvl_cn_sampleID), sampleID), setNames(all_of(gene_lvl_cn_geneID), geneID), everything())
      
      mut_mat <- bind_rows(mut_mat, gene_level_deletion_df)

    }
    
    
    # Format into a matrix
    collapse_string <- function(x){ 
      str_c(x, collapse = ",")
    }
    
    oncoprint_mat <- mut_mat %>%
      pivot_wider(names_from = !!sampleID, values_from = variant) %>% 
      rowwise() %>% mutate_all(collapse_string) 
    mat <- as.matrix(oncoprint_mat[,-1])
    rownames(mat) <- oncoprint_mat[[geneID]]
    mat[mat == ""] <- NA
    return(mat)
}

annotate_no_call <- function(x, samples){
  to_add <-  samples[!samples %in% colnames(x)]
  if(length(to_add) >= 1){
    missing <- matrix(data = NA, nrow = nrow(x), ncol = length(to_add))
    colnames(missing) = to_add
    rownames(missing) = rownames(x)
    output <- cbind(missing,x)
    output <- output[,samples]
  } else {
    output <- x[,samples]
  }
  return(output)
}
  