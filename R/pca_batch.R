#' @title PCA Plot with batch information
#'
#' @description Plots the first two principal components for all the samples of
#'  a gene expression dataset along with the batch information
#'
#' @param exprData A matrix containing gene expression data. Row names should be
#' samples and column names should be genes.
#' @param batch.info A data frame containing the samples names and details of the
#' batch they belong to.
#' @param batch Title of the batch the data is being corrected for.
#' @param NameString  string that will be appear in all output filenames. Default="" (empty string)
#' @param when String indicating for which dataset is the PCA plot being created
#'
#'  @import grDevices
#'  @import stats
#'  @import graphics
#'
#' @export
pca_batch <- function(exprData, batch.info, batch, NameString = "", when = ""){

  print("===========================Plotting PCs along with batch=====================")

  #matching sample IDs with batches
  if (is.character(batch.info[,1])){
    id <- as.character(rownames(exprData))
    s<- match (id, batch.info[,1])

  } else if (is.numeric(batch.info[,1])){
    id <-as.numeric(rownames(exprData))
    s<- match (id, batch.info[,1])
  }


  #calculating the principal components and adding the data to pca_data
  print("Calculating Principal Components...")
  pca.dat <- prcomp(exprData, center = TRUE, scale. = TRUE)
  print("PCs calculated")


  pca_data <- cbind.data.frame(pca.dat$x[, 1:2], batch.info[s,2])
  colnames(pca_data)[3] <- "Batch"
  pca_data[,3] <- as.factor(pca_data[,3])


  #plot PCA
  date <- as.character(format(Sys.Date(), "%Y%m%d"))
  plotFile <- ifelse (NameString =="",
                      paste0(date, "_plot_", batch, "_pca_", when, ".pdf"),
                      paste0(date, "_plot_", NameString, "_", batch, "_pca_", when, ".pdf"))
  #pdf (plotFile)
  pca_plot <- ggplot2::ggplot(data = pca_data) +

    ggplot2::geom_point(size =6,
                        alpha= 0.6,
                        ggplot2::aes(x=PC1, y=PC2, colour=Batch))+

    ggplot2::labs(title=paste("PCA with", batch, "information", sep = " "),
                  x = "PC1",
                  y = "PC2",
                  colour = "Batch")+

    ggplot2::theme(

      #Adjusting axis titles, lines and text
      axis.title = ggplot2::element_text(size = 15),
      axis.line = ggplot2::element_line(size=0.75),
      axis.text = ggplot2::element_text(size=15, colour ="black"),
      #Center align the title
      plot.title = ggplot2::element_text(face = "bold", hjust =0.5, size =20),

      #Adjust legend title and text
      legend.title = ggplot2::element_text(size = 15, face = "bold"),
      legend.text = ggplot2::element_text(size = 15),

      # Adjust panel border
      panel.border = ggplot2::element_rect(fill=NA, size= 0.75),

      # Remove panel grid lines
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),

      # Remove panel background
      panel.background = ggplot2::element_blank())+

     ggplot2::scale_colour_manual(values=c("#fcba03",  "#19e01c", "#ff470f",
                                           "#0fdbca", "#ff217e","#405ce6",
                                           "#6b6769","#b264ed"))
  return(pca_plot)

  #plot(pca_plot)
  #dev.off()

  #print(paste0("Plotted PC1 & PC2 with batch to ", plotFile))

}
