#' @title PCA plot with NMF cluster information
#'
#' @description Clusters the data using NMF (Non-negative matrix factorization)
#' after finding the optimal number of clusters in a dataset using the value of
#' cophenetic coefficients.The results of the clustering are used along with
#' PCA to see whether all the samples of a batch lie in the same cluster.
#' The cophenetic coefficient plot, the PCA biplot, and the files containing the
#' cophentic coefficeients (for number of clusters: 2 to 7) and the clustering
#' information (for optimal k) is saved to the NMF folder created in the working
#' directory.
#'
#'
#' @param expr gene expression dataset (rows should be genes, column should be samples)
#' @param batch.info contains the samples names and the batches they belong to
#' @param nrun The number of runs for NMF, the default number is 30. The consensus matrices from the NMF results are used to compute the cophenetic coefficients
#' @param batch title of the batch being used for correction
#' @param NameString  string that will be appear in all output filenames. Default="" (empty string)
#' @param when String indicating when the clustering is taking place (before batch correction or after batch correction?)
#'
#' @return Returns the optimal number of clusters (k) that has the maximum
#' average silhouette width (ignoring k=2) and was used for clustering and plotting.
#'
#' @import grDevices
#' @import utils
#' @import graphics
#'
#' @export
NMF_PCA <- function(expr, batch.info, nrun = 30, batch= "Batch", NameString = "", when){

  print(paste0("=========================== Performing NMF ", when, "====================="))
  #creating folder to store NMF data
  dir <- getwd()
  dir.create(paste0(dir, "/", "NMF_", when))

  #median centre data
  med <- apply (expr, MARGIN =1, median)
  median_centered_data<- expr - med

  #removing negative values for running NMF
  eps <- .Machine$double.eps
  median_centered_data[median_centered_data<=0 ] = eps
  median_centered_data <- as.matrix(median_centered_data)

  #running nmf
  NMF_estimate <- lapply(c(2:7), function(x){
    estimate <- NMF::nmf(median_centered_data, x, nrun = nrun, method = "lee", seed = 123)
    cluster <- cbind(names(NMF::predict(estimate)), NMF::predict(estimate))
    result <- list(estimate, cluster)
  })

  #unlist consensus matrices and clusters
  n <- dim(expr)[2]
  consensusMatrices=lapply(c(1:6),
                           function(x){s <- matrix(unlist(NMF_estimate[[x]][[1]]@consensus),
                                                   nrow =n, ncol = n)})
  #cophenetic coefficient
  coph <- data.frame(cluster=c(2:7), coph.coef = sapply(consensusMatrices, NMF::cophcor))

  #writing cophenetic coefficients to file
  date <- as.character(format(Sys.Date(), "%Y%m%d"))
  write.table(coph, file = ifelse(NameString == "",
                                  paste0(dir, "/", "NMF_", when, "/", date, "_", batch, "_", when,  "_cophenetic_coefficients.txt"),
                                  paste0(dir, "/", "NMF_", when, "/", date, "_", batch, "_", NameString, "_", when,  "_cophenetic_coefficients.txt")),
              sep = "\t")

  #finding the number of clusters with highest cophenetic coefficient
  num.opt <- coph[which(coph[,2]==max(coph[-1,2])),1]
  print(paste0("Number of clusters=", num.opt, " with the highest cophenetic coefficient = ", coph[which(coph[,1]==num.opt),2] ))

  #writing the clustering information to file
  cluster.info<-NMF_estimate[[num.opt-1]][[2]]
  colnames(cluster.info)<- c("Samples", "Class")
  write.table(cluster.info,
              file= ifelse(NameString == "",
                           paste0(dir, "/", "NMF_", when, "/", date, "_", when,  "_NMF_cluster_info.txt"),
                           paste0(dir, "/", "NMF_", when, "/", date, "_", NameString, "_", when,  "_NMF_cluster_info.txt")),
              sep = "\t", row.names = FALSE)

  #Calculating Principal Components
  print("Calculating Principal Components...")
  exprData <- t(expr)
  pca.data <- prcomp(exprData, center = TRUE, scale. = TRUE)

  #matching sample IDs with batches
  if (is.character(batch.info[,1])){
    id <- as.character(rownames(exprData))
    mat1<- match (id, batch.info[,1])
    mat2<- match (id, cluster.info[,1])

  } else if (is.numeric(batch.info[,1])){
    id <-as.numeric(rownames(exprData))
    mat1<- match (id, batch.info[,1])
    mat2<- match (id, as.numeric(cluster.info[,1]))
  } else{
    stop ("Sample IDs are neither characters nor numeric; Convert the sample IDs to either character or integers")
  }

  #combining PCs, batch information, K-means cluster information
  pca_data <- cbind.data.frame(pca.data$x[, 1:2], batch.info[mat1, 2], cluster.info[mat2, 2])
  colnames(pca_data)[3] <- "Batch"
  colnames(pca_data)[4] <- "NMF_cluster"

  pca_data[,3] <- as.factor(pca_data[,3])
  pca_data[,4] <- as.factor(pca_data[,4])

  #plotting cophenetic coefficient plot and PCA plot with batches and NMF information for the number of clusters with highest cophenetic coefficient
  if (NameString==""){
    plotFile <- paste0(date, "_plot_", batch, "_pca_with_NMF_info_", when, ".pdf")
  } else{
    plotFile <- paste0(date, "_plot_", NameString, "_", batch, "_pca_with_NMF_info_", when, ".pdf")
  }
  print(paste0("Plotting Cophenetic coefficient Plot and Principal Component Analysis biplot (with batches and clustering information) to the file ", plotFile))
  pdf (plotFile)

  #cophenetic coefficient plot
  coph.plot <- ggplot2::ggplot(data = coph, ggplot2::aes(x =cluster, y =coph.coef)) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = num.opt), size =1,
                        linetype = 'dashed', colour = "dodgerblue3")+
    ggplot2::geom_point(size = 2.5, colour = "orange")+
    ggplot2::geom_path(size = 1.25, colour = "orange") +
    ggplot2::labs(x = "Number of clusters",
                  y = "Cophenetic Coefficient",
                  title = "Optimal number of clusters for NMF clustering") +

    ggplot2::theme_classic()

  plot(coph.plot)

  #PCA with batch and kmeans info plot
  NMF_pca_plot <- ggplot2::ggplot(data = pca_data) +

    ggplot2::geom_point(ggplot2::aes(x=PC1,
                                     y=PC2,
                                     colour =Batch,
                                     shape =NMF_cluster), size=6, alpha = 0.6)+
    ggplot2::labs(title= paste0("PCA plot for ", batch, "; num of clusters= ", num.opt),
                  colour = "Batch",
                  shape = "NMF Cluster")+
    ggplot2::theme(

      #Adjusting axis titles, lines and text
      axis.title = ggplot2::element_text(size = 15),
      axis.line = ggplot2::element_line(size =0.75, colour = "black"),
      axis.text = ggplot2::element_text(size=15, colour ="black"),

      #Center align the title
      plot.title = ggplot2::element_text(face = "bold", hjust =0.5, size =20),

      #Adjust legend title and text
      legend.title = ggplot2::element_text(size = 10, face = "bold"),
      legend.text = ggplot2::element_text(size = 10),

      # Remove panel border
      panel.border = ggplot2::element_rect(fill=NA, size= 0.75),

      # Remove panel grid lines
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),

      # Remove panel background
      panel.background = ggplot2::element_blank()) +

    ggplot2::scale_color_manual(values=c("#fcba03",  "#19e01c", "#ff470f",
                                         "#0fdbca", "#ff217e","#405ce6",
                                         "#6b6769","#b264ed"))

  plot(NMF_pca_plot)

  dev.off()


}
