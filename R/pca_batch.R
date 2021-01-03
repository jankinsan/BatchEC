#' @title PCA Plot
#'
#' @description Plots the first two principal components for all the samples of
#'  a gene expression dataset along with the batch information
#'
#' @param exprData A matrix containing gene expression data. Row names should be
#' samples and column names should be genes.
#' @param batch.info A data frame containing the samples names and details of the
#' batch they belong to.
#' @param batch Title of the batch the data is being corrected for.
#' @param plotFile Name of the file in which the PCA is to be plotted.
#'
#'  @import grDevices
#'  @import stats
#'  @import graphics
#'
#' @export
pca_batch <- function(exprData, batch.info, batch, plotFile){

  print("===========================Plotting PCs along with batch=====================")

  #matching sample IDs with batches
  id <- as.numeric(rownames(exprData))
  s<- match (id, batch.info[,1])

  #calculating the principal components and adding the data to pca_data
  print("Calculating Principal Components...")
  pca1 <- prcomp(exprData, center = TRUE, scale. = TRUE)


  pca_data <- cbind.data.frame(pca1$x[, 1:2], batch.info[s,2])
  colnames(pca_data)[3] <- "Batch"
  pca_data[,3] <- as.factor(pca_data[,3])


  #plot PCA
  date <- as.character(format(Sys.Date(), "%Y%m%d"))
  plotFile <- paste0(date, "_", plotFile)
  pdf (plotFile)
  pca_plot <- ggplot2::ggplot(data = pca_data) +

    ggplot2::geom_point(ggplot2::aes(x=PC1, y=PC2, colour =Batch))+

    ggplot2::labs(title=paste("PCA with batch for", batch, sep = " "),
                  x = "PC1",
                  y = "PC2",
                  colour = batch)

    ggplot2::theme(

      #adjusting axis lines and text
      axis.line.x = ggplot2::element_line(size =0.75),
      axis.line.y = ggplot2::element_line(size =0.75),
      axis.text.x = ggplot2::element_text(angle = 90, size=8.5, colour ="black"),
      axis.text.y = ggplot2::element_text(size=8.5, colour ="black"),

      #Center align the title
      plot.title = ggplot2::element_text(face = "bold", hjust =0.5),

      #remove legend title
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 10),

      # Adjust panel border
      panel.border = ggplot2::element_blank(),

      # Remove panel grid lines
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),

      # Remove panel background
      panel.background = ggplot2::element_blank()


  )

  plot(pca_plot)
  dev.off()

  print(paste0("Plotted PC1 & PC2 with batch in ", plotFile))

}

