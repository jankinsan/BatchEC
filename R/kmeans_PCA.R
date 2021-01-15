#' @title PCA plot with k-means information
#'
#' @description Clusters the data using k-means after finding the optimum number
#' of clusters using the silhouette method and displays the results on a PCA plot
#' saved to a file. The files containing the avg. silhouette width (for k = 2 to
#' k = 7) and the cluster information re saved to the kmeans folder within the working directory
#'
#'
#' @param exprData gene expression dataset (rows should be samples, column should be genes)
#' @param batch.info contains the samples names and the batches they belong to
#' @param batch title of the batch being used for correction
#' @param NameString  string that will be appear in all output filenames. Default="" (empty string)
#' @param when String indicating when the clustering is taking place (before batch correction or after batch correction?)
#'
#' @return Returns the optimal number of clusters (k) that has the maximum average silhouette width
#' and was used for clustering
#'
#' @import grDevices
#' @import stats
#' @import utils
#' @import graphics
#'
#' @export
kmeans_PCA <- function(exprData, batch.info, batch= "Batch", NameString = "", when){

  print(paste0("=========================== Performing k-means ", when, "====================="))

  #calculating distance matrix
  print ("Calculating distance matrix")
  distMatrix <- dist(exprData, method = "euclidean")

  #using the silhouette method for determining the optimal number of clusters for k-means
  print ("Determining the optimal number of clusters for k-means...")
  silh <- data.frame()
  for (k in 2:7){
  cluster<- kmeans(exprData, centers = k, iter.max = 100, nstart = 25,
                           set.seed(12345))
  sil.data <- cluster::silhouette(cluster$cluster, distMatrix)
  value <- c(k, mean(sil.data[,3]))
  silh <- rbind (silh, value)
  }

  colnames(silh)<- c("k", "silh.width")

  #writing avg sil width to file and plotting the values
  date <- as.character(format(Sys.Date(), "%Y%m%d"))
  write.table(silh, file = paste0(date, "_", NameString, "_", when, "_avg_silhouette_width_k-means.txt"),
              sep = "\t")

  plot_silh_file <- paste0(date, "_", NameString, "_", when, "_plot_optimal_clusters_silhouette_method.pdf")
  pdf (plot_silh_file)
  print(factoextra::fviz_nbclust(exprData, kmeans, method = "silhouette"))
  dev.off(which=dev.cur())

  #finding the k with maximum avg. silhouette width
  max.sil <- silh[which(silh[,2]==max(silh[,2])),]
  print(paste0("k=", max.sil[1], " has the maximum Average Silhouette Width = ",
               max.sil[2]))
  k <- as.numeric(max.sil[1])


  #clustering for optimal k
  print("Clustering data ...")
  clustered_data <- kmeans(exprData, centers = k, iter.max = 100, nstart = 25,
                           set.seed(12345))
  cluster.info <- cbind(names(clustered_data$cluster), clustered_data$cluster)
  colnames(cluster.info)<- c("Samples", "Cluster")
  write.table(cluster.info,
              file= paste0(date, "_", NameString, "_", when,  "_k-means_cluster_info.txt"),
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
  }



  #combining PCs, batch information, K-means cluster information
  pca_data <- cbind.data.frame(pca.data$x[, 1:2], batch.info[mat1, 2], cluster.info[mat2, 2])
  colnames(pca_data)[3] <- "Batch"
  colnames(pca_data)[4] <- "kmeans_cluster"

  pca_data[,3] <- as.factor(pca_data[,3])
  pca_data[,4] <- as.factor(pca_data[,4])

  #plotting PCA with batch along with kmeans information for k
  plotFile <- paste0(date, "_", NameString, "_", batch, "_pca_with_kmeans_info_", when, ".pdf")

  print(paste0("Plotting Principal Components along with clustering information to the file ", plotFile))
  pdf (plotFile)
  kmeans_pca_plot <- ggplot2::ggplot(data = pca_data) +

    ggplot2::geom_point(ggplot2::aes(x=PC1,
                   y=PC2,
                   colour =Batch,
                   shape =kmeans_cluster), size=2)+
    ggplot2::labs(title= paste0("PCA plot for ", batch, "; k = ", k),
                  colour = "Batch",
                  shape = "k-means Cluster")+
    ggplot2::theme(

      #adjusting axis lines and text
      axis.line.x = ggplot2::element_line(size =0.75),
      axis.line.y = ggplot2::element_line(size =0.75),
      axis.text.x = ggplot2::element_text(angle = 90, size=8.5, colour ="black"),
      axis.text.y = ggplot2::element_text(size=8.5, colour ="black"),

      #Center align the title
      plot.title = ggplot2::element_text(face = "bold", hjust =0.5),

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

  plot(kmeans_pca_plot)

  dev.off()

  return (k)
}

