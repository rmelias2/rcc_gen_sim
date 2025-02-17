#' @param plot Dataframe containing plot information
#' @param x X-axis variable
#' @param y Y-axis variable
#' @param color_var Variable to color by
#' @param shape_var Variable to shape by
#' @param labels Genes to label
#' @param labels_value Column to use for gene labels

ggplot_alteration_frequency_scatter <- function(plot, x_val, y_val, color_var, shape_var = NULL, labels = NULL, xlim = 1, ylim = 1, labels_value = "gene", ...) {
  required_packages <- c("ggplot2", "dplyr", "ggrepel", "ggpubr", "RColorBrewer")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("Package", pkg, "is required but not installed. Please install it using install.packages()."))
    } else {
      library(pkg, character.only = TRUE)
    }
  }
  
  g <- ggplot(plot, aes(x = .data[[x_val]], y = .data[[y_val]])) +   geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed", linewidth = 0.3) 
  
  if (!is.null(labels)) {
    if (is.null(shape_var)) {
      g <- g + 
        geom_jitter(data = plot[!(plot[[labels_value]] %in% labels), ], aes(color = .data[[color_var]]), width = 0.02, height = 0.02) +
        geom_point(data = plot[plot[[labels_value]] %in% labels, ], aes(color = .data[[color_var]])) +
        geom_label_repel(data = plot[plot[[labels_value]] %in% labels, ],
                         aes(label = .data[[labels_value]]),     
                         color = "black", 
                         force = 5, 
                         segment.size = 0.5,
                         segment.curvature = 0,
                         box.padding = 1, 
                         nudge_x = 0.2,
                         nudge_y = 0.1,
                         point.size = 0.01,
                         min.segment.length = 0.5, 
                         max.overlaps = 100, size = 3)
    } else {
      g <- g + 
        geom_jitter(data = plot[!(plot[[labels_value]] %in% labels), ], aes(color = .data[[color_var]], shape = .data[[shape_var]]), width = 0.02, height = 0.02) +
        geom_point(data = plot[plot[[labels_value]] %in% labels, ], aes(color = .data[[color_var]], shape = .data[[shape_var]])) +
        geom_label_repel(data = plot[plot[[labels_value]] %in% labels, ],
                         aes(label = .data[[labels_value]]),     
                         color = "black", 
                         force = 5, 
                         segment.size = 0.5,
                         segment.curvature = 0,
                         box.padding = 1, 
                         nudge_x = 0.2,
                         nudge_y = 0.1,
                         point.size = 0.01,
                         min.segment.length = 0.5, 
                         max.overlaps = 100, size = 3)
    }
  } else {
    if (is.null(shape_var)) {
      g <- g + geom_jitter(aes(color = .data[[color_var]]), width = 0.02, height = 0.02)
    } else {
      g <- g + geom_jitter(aes(color = .data[[color_var]], shape = .data[[shape_var]]), width = 0.02, height = 0.02)
    }
  }
  
  # Assign Color Scale
  if (class(plot[[color_var]]) %in% c("character", "factor")) {
    unique_values <- unique(plot[[color_var]])
    palette <- brewer.pal(n = min(length(unique_values) - 1, 9), name = "Set1")
    colors <- setNames(c(palette, "gray"), c(unique_values[unique_values != "other"], "other"))
    
    g <- g + scale_color_manual(values = colors) 
  } else {
    myPalette <- colorRampPalette(brewer.pal(11, "RdYlBu"))
    #g <- g + scale_colour_steps2(low = "red", mid = "gray", high = "black", limits = c(0, 1), midpoint = 0.5)
    #g <- g + scale_colour_gradient(limits = c(0,1), low = "red", high = "black", breaks = c(0,0.1, 1))
    g <- g + scale_colour_gradientn(colours = myPalette(100), values=seq(0, 1, length.out=20), limits = c(0,1))
  }
  
  g <- g + 
  
    scale_x_continuous(limits = c(0, as.numeric(xlim)), oob = scales::squish) + 
    scale_y_continuous(limits = c(0, as.numeric(ylim)), oob = scales::squish) + 
    theme_minimal(base_size = 12) + 
    labs(...) 
  
  return(g)
}
