# HEATMAP/ONCOPRINT HELPER FUNCTIONS
library(circlize)

#' ## Stored Top Annotations for Oncoprints and Heatmaps
AFR.col = colorRamp2(c(0, 0.3, 0.5, 0.7, 1), 
                     c("#ffffb2", "#fecc5c", "#fd8d3c", "#f03b20", "#bd0026"))

EUR.col = colorRamp2(c(0, 0.3, 0.5, 0.7, 1), 
                     c("#eff3ff", "#bdd7e7", "#6baed6", "#3182bd", "#08519c"))

Error.col = colorRamp2(c(38, 60, 95), 
                      c("darkgreen", "#bae4bc", "#fd8d3c"))

lib_size.col = colorRamp2(c(1238517,5000000, 10000000 ,14054983), c("#e0f3db","#bae4bc", "#7bccc4", "#2b8cbe"))


Type.col <- c("Benign" = "green3", 
              "Tumor" = "red3")

pT.col = c("N" = "lightgray",
                "pT1" = "yellow3",
                "pT3" = "red4")



heatmap_legend_param_rna = list(title = "Normalized Expression",
                                title_gp = gpar(fontsize = 4, fontface = "bold"), grid_height = unit(2, "mm"),
                                grid_width = unit(2, "mm"), labels_gp = gpar(fontsize = 4))


annotation_legend_param = list(title_gp = gpar(fontsize = 4, fontface = "bold"), grid_height = unit(2, "mm"),
                               grid_width = unit(2, "mm"), labels_gp = gpar(fontsize = 4), annotation_name_align = TRUE)


heatmap_legend_param_cna = list(title = "Log2Fold Change",
                                title_gp = gpar(fontsize = 4, fontface = "bold"), grid_height = unit(2, "mm"),
                                grid_width = unit(2, "mm"), labels_gp = gpar(fontsize = 4))
