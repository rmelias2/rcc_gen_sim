#' Function dea_data_prep
#' Take summarized experiment object, subset based on variables of interest, and remove low count genes. 
#' @param se.object Summarized experiment object with an assay named `counts`
#' @param filter.var vector of variable names to filter
#' @param filter.values List of values for each variable in filter var on which to subset the data
#' @param group.var Variable name on which to assign groups for EdgeR low cound removal
#' 
#' ## Note on  `filterByExpr()`:
#' This function implements the filtering strategy that was described informally by Chen et al (2016).
#' Roughly speaking, the strategy keeps genes that have at least min.count reads in a worthwhile number samples.
#' More precisely, the filtering keeps genes that have CPM >= CPM.cutoff in MinSampleSize samples, where CPM.cutoff = min.count/median(lib.size)*1e6 and
#' MinSampleSize is the smallest group sample size or, more generally, the minimum inverse leverage computed from the design matrix.
#' If all the group samples sizes are large, then the above filtering rule is relaxed slightly. If MinSampleSize > large.n, then genes are kept if CPM >= CPM.cutoff in k samples where k = large.n + (MinSampleSize - large.n) * min.prop. This rule requires that genes are expressed in at least min.prop * MinSampleSize samples, even when MinSampleSize is large.
#' In addition, each kept gene is required to have at least min.total.count reads across all the samples. 
#' 
#' Default Params" `min.count = 10, min.total.count = 15, large.n = 10, min.prop = 0.7`
#' 
#' By adjusting the parameters above (and since all group sizes are larger than n, this keeps all genes with a minumum
#' CPM.cutoff in at least 1 sample. (2% of the smallest GGG group (intermediate, n = 31) < 1))



dea_data_prep <- function(se.object, filter.var = NULL, filter.values = NULL, group.var){
  required_packages <- c("DESeq2", "edgeR")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("Package", pkg, "is required but not installed. Please install it using install.packages()."))
    } else {
      library(pkg, character.only = TRUE)
    }
  }
    
  meta <- as.data.frame(colData(se.object))
  
  if(!is.null(filter.var)){
    for(variable in 1:length(filter.var)){
      meta <- meta %>% filter(.data[[sym(filter.var[variable])]] %in% filter.values[[variable]])
    }
  }
  # Subset for samples in the meta data file. Note, this code is dependent on a variable called "sampleID"
  se <- se.object[,se.object$sampleID %in% meta$sampleID]
  
  #Subset low count genes. Note, also dependent on an assay called "counts"
  dge <- DGEList(assay(se, "counts"), group = meta[[group.var]])
  keep <- filterByExpr(dge, min.count = 5, min.total.count = 15, large.n = 10, min.prop = 0.02)
  se <- se[keep,]
  return(se)
  
  }
