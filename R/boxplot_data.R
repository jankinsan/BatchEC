#' @title Boxplots showing Gene Expression before and after correction
#'
#' @description Plots boxplots showing the change in gene expression after
#' adjustment for batch effects using ComBat
#'
#' @param exp1 A dataframe containing gene expression values before batch correction.
#' Rows are genes, columns are samples.
#' @param exp2 A dataframe containing gene expression values after adjustment
#' for batch effects. #' Rows are genes, columns are samples.
#' @param NameString  String that will be added in all output filenames. Default=NA.
#' @param batch The title of the batch for which you want to evaluate and do the
#'  correction. Default = "Batch"
#'
#'  @import grDevices
#'
#'  @export
boxplot_data <- function(expr1, expr2, NameString, batch = "Batch"){

  print ("===========================Plotting boxplots showing Gene Expression=====================")

  #plotting gene expression before batch correction--------------------------
  data1 <- reshape2::melt(expr1)
  boxplot_before_correction <- ggplot2::ggplot(data1, ggplot2::aes(x =data1[["variable"]],
                                                                   y =data1[["value"]])) +
    ggplot2::geom_boxplot(color = "blue4",
                          inherit.aes = TRUE,
                          outlier.size = 0.5,
                          outlier.color = "limegreen")+
    ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(check.overlap = TRUE))+

    ggplot2::labs(title = "Gene Expression before batch correction",
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

  #plotting gene expression after batch correction--------------------------
  data2 <- reshape2::melt(expr2)
  boxplot_after_correction <- ggplot2::ggplot(data2, ggplot2::aes(x=data2[["variable"]],
                                                                   y=data2[["value"]])) +
    ggplot2::geom_boxplot(color = "red",
                          inherit.aes = TRUE,
                          outlier.size = 0.5,
                          outlier.color = "limegreen")+
    ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(check.overlap = TRUE))+

    ggplot2::labs(title = "Gene Expression after batch correction",
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


  date <- as.character(format(Sys.Date(), "%Y%m%d"))
  if (NameString == ""){
    plotFile <- paste0(date, "_", batch, "_gene_expression_boxplots.pdf")
  }else {
    plotFile <- paste0(date, "_", NameString, "_", batch, "_gene_expression_boxplots.pdf")
  }

  #plotting both the boxplots in a pdf file
  pdf(plotFile)
  plot (boxplot_before_correction)
  plot (boxplot_after_correction)
  print("Boxplots Plotted...")
  dev.off()
}
