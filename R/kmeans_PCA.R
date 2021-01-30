#' @title PCA plot with k-means information
#'
#' @description Clusters the data using k-means clustering after finding the
#' optimal number of clusters for the dataset using the silhouette method. The
#' results of the clustering are used along with PCA to see whether all the
#' samples of a batch lie in the same cluster. The silhouette plot, the PCA
#' biplot, and the files containing the avg. silhouette width (for k = 2 to
#' k = 7) and the clustering information (for optimal k) is saved to the k-means
#' folder created in the working directory.
#'
#'
#' @param exprData gene expression dataset (rows should be samples, column should be genes)
#' @param batch.info contains the samples names and the batches they belong to
#' @param batch title of the batch being used for correction
#' @param NameString  string that will be appear in all output filenames. Default="" (empty string)
#' @param when String indicating when the clustering is taking place (before batch correction or after batch correction?)
#'
#' @return Returns the optimal number of clusters (k) that has the maximum
#' average silhouette width (ignoring k=2) and was used for clustering and plotting.
#'
#' @import grDevices
#' @import stats
#' @import utils
#' @import graphics
#'
#' @export
kmeans_PCA <- function(exprData, batch.info, batch= "Batch", NameString = "", when){

  print(paste0("=========================== Performing k-means ", when, "====================="))
  #creating folder to store k-means data
  dir <- getwd()
  dir.create(paste0(dir, "/", "kmeans_", when))

  #calculating distance matrix
  print ("Calculating distance matrix")
  distMatrix <- dist(exprData, method = "euclidean")

  #using the silhouette method for determining the optimal number of clusters for k-means
  print ("Determining the optimal number of clusters for k-means and clustering data...")
  #clustering data for various k
  cluster_data<- lapply(c(2:7), function(x){kmeans(exprData,
                                                   centers= x,
                                                   iter.max = 1000,
                                                   nstart=25,
                                                   set.seed(12345))})


  silh <- as.data.frame(lapply(c(2:7),  function(x){
  c(x, mean(cluster::silhouette(cluster_data[[x-1]]$cluster, distMatrix)[,3]))}),
  col.names = c(2:7), row.names = c("k", "silWidth"))

  silh <- as.data.frame(t(silh))

  #writing avg sil width to file
  date <- as.character(format(Sys.Date(), "%Y%m%d"))
  write.table(t(silh), file = paste0(date, "_", NameString, "_", when, "_avg_silhouette_width_k-means.txt"),
              sep = "\t")


  #finding the k with maximum avg. silhouette width
  max.sil <- silh[1+which(silh[2:dim(silh)[1],2]==max(silh[2:dim(silh)[1],2])),]
  print(paste0("k=", max.sil[1], " has Average Silhouette Width = ", max.sil[2]))
  opt.k <- as.numeric(max.sil[1])

  #writing the clustering information to file
  clustered_data<-cluster_data[[opt.k - 1]]
  cluster.info <- cbind(names(clustered_data$cluster), clustered_data$cluster)
  colnames(cluster.info)<- c("Samples", "Cluster")
  write.table(cluster.info,
              file= ifelse(NameString == "",
                           paste0(dir, "/", "kmeans_", when, "/", date, "_", when,  "_k-means_cluster_info.txt"),
                           paste0(dir, "/", "kmeans_", when, "/", date, "_", NameString, "_", when,  "_k-means_cluster_info.txt")),
              sep = "\t")

  #Calculating Principal Components
  print("Calculating Principal Components...")
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
  colnames(pca_data)[4] <- "kmeans_cluster"

  pca_data[,3] <- as.factor(pca_data[,3])
  pca_data[,4] <- as.factor(pca_data[,4])

  #plotting silhouette plot and PCA plot with batches and kmeans information for opt.k
  if (NameString==""){
    plotFile <- paste0(date, "_", batch, "_pca_with_kmeans_info_", when, ".pdf")
  } else{
    plotFile <- paste0(date, "_", NameString, "_", batch, "_pca_with_kmeans_info_", when, ".pdf")
  }
  print(paste0("Plotting Silhouette Plot and Principal Component Analysis biplot (with batches and clustering information) to the file ", plotFile))
  pdf (plotFile)

  #silhouette plot
  silh.plot <- ggplot2::ggplot(data = silh, ggplot2::aes(x =k, y =silWidth)) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = opt.k), size =1,
                        linetype = 'dashed', colour = "dodgerblue3")+
    ggplot2::geom_point(size = 2.5, colour = "orange")+
    ggplot2::geom_path(size = 1.25, colour = "orange") +
    ggplot2::labs(x = "Number of clusters, k",
                  y = "Average Silhouette Width",
                  title = "Optimal number of clusters for k-means") +

    ggplot2::theme_classic()

  plot(silh.plot)

  #PCA with batch and kmeans info plot
  kmeans_pca_plot <- ggplot2::ggplot(data = pca_data) +

    ggplot2::geom_point(ggplot2::aes(x=PC1,
                   y=PC2,
                   colour =Batch,
                   shape =kmeans_cluster), size=6, alpha = 0.6)+
    ggplot2::labs(title= paste0("PCA plot for ", batch, "; k = ", opt.k),
                  colour = "Batch",
                  shape = "k-means Cluster")+
    ggplot2::theme(

      #adjusting axis titles, lines and text
      axis.title = ggplot2::element_text(size = 15),
      axis.line = ggplot2::element_line(size =0.75, colour = "black"),
      axis.text = ggplot2::element_text(size=15, colour ="black"),

      #Center align the title
      plot.title = ggplot2::element_text(face = "bold", hjust =0.5, size =20),

      # Remove panel border
      panel.border = ggplot2::element_blank(),

      # Remove panel grid lines
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),

      # Remove panel background
      panel.background = ggplot2::element_blank()
    ) +
    scale_fill_manual(values=c("#000000", "#E69F00", "#56B4E9", "#009E73",
                               "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))

  plot(kmeans_pca_plot)

  dev.off()

  return (k)
}


