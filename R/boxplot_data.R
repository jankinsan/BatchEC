#' @title Plots boxplots showing Gene Expression across samples
#'
#' @description Plots boxplots that can be used to see the change in gene expression after
#' adjustment for batch effects using ComBat
#'
#' @param expr A dataframe containing gene expression values.
#' Rows are genes, columns are samples. Could be the dataframe before or after correction
#' @param NameString  String that will be added in all output filenames. Default=NA.
#' @param batch The title of the batch for which you want to evaluate and do the
#'  correction. Default = "Batch"
#' @param when String indicating if the expression data is from before or after correction. Use "before" or "after"
#'
#' @import grDevices
#'
#' @export
boxplot_data <- function(expr, when, NameString, batch = "Batch"){

  print(paste0("================Plotting boxplot ", when, "====================="))

  #plotting gene expression
  data <- reshape2::melt(expr)

  if(when == "before"){
  boxplot_expression <- ggplot2::ggplot(data, ggplot2::aes(x=variable,
                                                          y= value)) +
    ggplot2::geom_boxplot(color = "blue4",
                          inherit.aes = TRUE,
                          outlier.size = 0.5,
                          outlier.color = "limegreen")+
    ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(check.overlap = TRUE))+

    ggplot2::labs(title = paste0("Gene Expression ", when, " Correction for ", batch),
                  x = "Samples",
                  y = "Expression") +

    ggplot2::theme(
      #adjusting axis lines and text
      axis.line.x = ggplot2::element_line(size =0.5),
      axis.line.y = ggplot2::element_line(size =0.5),
      axis.text.x = ggplot2::element_text(size =10, colour ="black"),
      axis.text.y = ggplot2::element_text(size=10, colour ="black"),

      #Center align the title
      plot.title = ggplot2::element_text(face = "bold", hjust =0.5),

      # Remove panel grid lines
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),

      # Remove panel background
      panel.background = ggplot2::element_blank())
  }
  else if (when == "after"){
    boxplot_expression <- ggplot2::ggplot(data, ggplot2::aes(x =variable,
                                                              y =value)) +
      ggplot2::geom_boxplot(color = "red",
                            inherit.aes = TRUE,
                            outlier.size = 0.5,
                            outlier.color = "limegreen")+
      ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(check.overlap = TRUE))+

      ggplot2::labs(title = paste0("Gene Expression ", when, " after Batch Correction"),
                    x = "Samples",
                    y = "Expression") +

      ggplot2::theme(
        #adjusting axis lines and text
        axis.line.x = ggplot2::element_line(size =0.5),
        axis.line.y = ggplot2::element_line(size =0.5),
        axis.text.x = ggplot2::element_text(size =10, colour ="black"),
        axis.text.y = ggplot2::element_text(size=10, colour ="black"),

        #Center align the title
        plot.title = ggplot2::element_text(face = "bold", hjust =0.5),

        # Remove panel grid lines
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),

        # Remove panel background
        panel.background = ggplot2::element_blank())
  } else{
    stop("Please specify the string 'when', should be either 'before' or 'after'")
  }

  date <- as.character(format(Sys.Date(), "%Y%m%d"))
  if (NameString == ""){
    plotFile <- paste0(date, "_", batch, "_gene_expression_boxplot_", when,".pdf")
  }else {
    plotFile <- paste0(date, "_", NameString, "_", batch, "_gene_expression_boxplot_", when,".pdf")
  }

  #plotting both the boxplots in a pdf file
  pdf(plotFile, height= 6.5, width = 13)
  plot (boxplot_expression)
  dev.off()
  print(paste0("Boxplot ", when, " Batch Correction plotted to ", plotFile))
}
