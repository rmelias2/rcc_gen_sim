#' Function multi_model_dea
#' @param deseq2.output Summarized experiment object with an assay named `counts`
#' @param exp.name vector of variable names to filter
#' @param additional.comparisons List of additional comparisons. Takes the format of list("comparison_name" = c("variable", "comparison_value", "ref_value")).
#' @param save.default Default T. Will save the default comparisons of DESEQ2 less then intercept. If False, will only save the comparisons in `additional.comparisons` 
#' The list names are used to name the output. Default value is null. Default comparisons can be identified using `resultsNames(deseq2.output)[-1]`


dea_save_output <- function(deseq2.output, exp.name, additional.comparisons = NULL, save.default = T){
  dds = deseq2.output
  #Make sub directory for output
  output.directory <-  here(output.path, exp.name)
  if (!dir.exists(here(output.directory))) dir.create(here(output.directory))
  
  #Save RDS Object
  write_rds(dds, here(output.directory, "DEseq2_object.rds"))
  
  
  #Save .csv of outputs
  maxCooks <- apply(assays(dds)[["cooks"]], 1, max)
  
  #Save Defaults
  if(save.default){
  dds_names <- resultsNames(dds)[-1]
  for (i in seq_along(dds_names)){
    output <- as.data.frame(results(dds, name = dds_names[i]))
    output$Gene_symbol <- rownames(results(dds))
    output$max_cooks <- maxCooks
    write_tsv(output, here(output.directory, paste0("Results_", dds_names[i], ".txt")))
  }
  }
  ## Additional Comparisons: 
  if(!is.null(additional.comparisons)){
    ac = additional.comparisons
    for(i in seq_along(1:length(ac))){
      output <- as.data.frame(results(dds, contrast = ac[[i]]))
      output <- output %>% mutate(Gene_symbol = rownames(output),
                                  max_cooks = maxCooks)
      write_tsv(output, here(output.directory, paste0("Results_", names(ac)[i], ".txt")))
    }
  }
}