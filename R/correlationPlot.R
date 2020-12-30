#' @title Correlation Scatter Plot
#'
#' @description The function calculates the correlation coefficient using
#' the expression data before and after batch correction and also plots a scatter
#' plot to show the correlation.
#'
#' @param data1 A matrix containing the gene expression data before batch correction
#' @param data2 A matrix containing the gene expression data after batch correction
#' @param fileName String to be added to the name of the file to which correlation
#'  scatter plot would be saved
#' @param batch The tile of the batch
#'
#' @import grDevices
#'
#' @export
correlationPlot<- function(data1, data2, fileName="", batch){

  print("===========================Plotting scatter plot to show correlation=====================")

  #matching sample IDs and Genes
  mat_gene <- match(colnames(data1), colnames(data2))
  mat_row <- match(row.names(data1), row.names(data2))

  data2 <- data2[mat_row, mat_gene]

  before_correction<- c(as.numeric(data1))
  after_correction<- c(as.numeric(data2))

  #calculating correlation
  corr_coef <- cor(before_correction, after_correction)

  #plotting gene expression values before and after batch correction
  x = cbind.data.frame(before_correction, after_correction)

  date <- as.character(format(Sys.Date(), "%Y%m%d"))
  if  (fileName ==""){
    plotFile <- paste0(date, "_plot_", batch,"_correlationPlot.pdf")
  } else{
    plotFile <- paste0(date, "_plot_", fileName, "_", batch,"_correlationPlot.pdf")
  }
  print(paste0("Correlation scatter plot saved to: ", plotFile))
  pdf (plotFile)

  plot1 <- ggplot2::ggplot(x, ggplot2::aes(x=before_correction, y=after_correction))+
    ggplot2::geom_point() +
    ggplot2::labs(x = "Gene Expression before Correction",
         y = "Gene Expression after Correction",
         title = "Sactterplot showing Correlation between data before and after batch correction",
         subtitle = paste("Pearson Correlation coefficient", "=", corr_coef, sep= " "))+
    ggplot2::theme(

      #adjusting axis lines and text
      axis.line.x = ggplot2::element_line(size =0.75),
      axis.line.y = ggplot2::element_line(size =0.75),
      axis.text.x = ggplot2::element_text(angle = 90, size=8.5, colour ="black"),
      axis.text.y = ggplot2::element_text(size=8.5, colour ="black"),

      #Center align the title
      plot.title = ggplot2::element_text(face = "bold"),

      # Remove panel border
      panel.border = ggplot2::element_blank(),

      # Remove panel grid lines
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),

      # Remove panel background
      panel.background = ggplot2::element_blank(),

      # Add axis line
      axis.line = ggplot2::element_line(colour = "grey")

    )


  plot(plot1)
  dev.off()

}
