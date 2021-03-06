#' @title Batch Effects Evaluation and Correction
#'
#'
#' @description Evaluates if the batch effects with a gene expression dataset using linear
#' regression and if the batch is associated with the data, the data is batch-corrected
#' using the ComBat Algorithm.
#'
#' Please set your working directory before you call the function as all output files will be saved to this folder.
#'
#'
#'
#' @param expr1 gene expression data; rows should be genes and columns should be samples.
#' @param batch.info contains the batch information. The first column should be sample names.
#' @param batch The title of the batch for which you want to evaluate and do the
#'  correction. Default = "Batch"
#' @param NameString  String that will be added in all output filenames. Default=NA.
#' @param discrete.batch Logical value indicating whether the samples are already
#' grouped in discrete batches. If the value is FALSE, contiguous batch information is
#' clustered into discrete batches. Useful for clustering batch variables that are
#' contiguous, for example, used reads and useful reads in mClust. Default = TRUE
#' @param clus.method Method to be used for clustering. "km" denotes k-means clustering, "NMF"
#' denotes NMF with 30 runs. Default employs both the methods, using *clus.method= c("NMF", "km")*.
#' @param nrun.NMF number of runs for NMF. Default = 30.
#'
#' @import utils
#' @import mclust
#'
#' @examples
#' \dontrun{
#' beacon(expr1 = dataset,
#'        batch.info = batch_info,
#'        batch = "Batch",
#'        discrete.batch = TRUE)
#' }
#'
#'
#' @export
beacon <- function(expr1, batch.info, batch, NameString = "", discrete.batch = TRUE, clus.method = c("NMF", "km"), nrun.NMF = 30){

  #matching batch
  mat <- match(batch, colnames(batch.info))
  batch.info <- batch.info [,c(1, mat)]

  if (discrete.batch== FALSE){
    batch.info[,2] <- mclust_cluster(batch.info[,2])

    #writing the batch cluster information to file
    date <- as.character(format(Sys.Date(), "%Y%m%d"))
    clusterInfoFile <- paste0(date, "_", NameString, "_", batch, "_batch_info_mClust.txt")
    print(paste0("Writing the batch information from mClust for ", batch, " to: ", clusterInfoFile))
    write.table(batch.info,
                file = clusterInfoFile,
                sep = "\t",
                col.names = TRUE,
                row.names = FALSE)
  }
  batch.info[,2]<- as.factor(batch.info[,2])


  #removing genes with zero variance
  std_genes <- apply(expr1, MARGIN =1, sd)
  genes_sd_0 <- length(which (std_genes ==0))
  remgenes <- length (std_genes) - genes_sd_0
  expr1 <- expr1[which (std_genes > 0),]
  print(paste0("Removed ", genes_sd_0, " genes with zero variance..."))
  print (paste0(remgenes, " genes remain..."))


  exprData1 <- t(expr1)

  #linear regression before batch correction
  p_val_before <- lin_reg(exprData=exprData1,
                          batch.info = batch.info,
                          batch=batch,
                          NameString = NameString,
                          when = "before")

  #checking if any of the p values are less than 0.05
  if(all(p_val_before[,2]>0.05)){
    stop("Halting execution since the batch specified is not associated with data...")

  }else{

    print("Batch is associated with the data...")

    #pca with batch before correction
    pca_batch (exprData = exprData1,
               batch.info= batch.info,
               batch= batch,
               NameString = NameString,
               when = "before_ComBat_correction")

    if(any(clus.method=="km")){
    #pca with batch and kmeans before batch correction
      k_before <- kmeans_PCA(exprData = exprData1,
                             batch.info = batch.info,
                             batch = batch,
                             NameString = NameString,
                             when = "before_correction")
    }
    if(any(clus.method=="NMF")){
      #NMF with PCA before correction
      NMF_PCA(expr=expr1,
              batch.info = batch.info,
              nrun = nrun.NMF,
              batch = batch,
              NameString = NameString,
              when = "before_correction")
    }

    #Batch Correction using ComBat
    expr2 <- ComBat_data (expr = expr1,
                          batch.info= batch.info,
                          batch = batch,
                          NameString = NameString)
    exprData2 <- t(expr2)

    #linear regression after batch correction
    p_val_after <- lin_reg(exprData=exprData2,
                           batch.info = batch.info,
                           batch=batch,
                           NameString = NameString,
                           when = "after")

    #pca with batch after correction
    pca_batch (exprData = exprData2,
               batch.info= batch.info,
               batch= batch,
               NameString = NameString,
               when = "after_ComBat_correction")

    if(any(clus.method=="km")){
    #pca with batch and k-means after correction
    k_after <- kmeans_PCA(exprData = exprData2,
                          batch.info = batch.info,
                          batch = batch,
                          NameString = NameString,
                          when = "after_correction")
    }
    if(any(clus.method=="NMF")){
      #NMF with PCA after correction
      NMF_PCA(expr=expr2,
              batch.info = batch.info,
              nrun = nrun.NMF,
              batch = batch,
              NameString = NameString,
              when = "after_correction")
    }

    #PCA Proportion of Variation
    pca_prop_var(batch.title = batch,
                 plotFile = ifelse(NameString == "",
                                   paste0("plot_pca_prop_var_", batch, ".pdf"),
                                   paste0(NameString, "_plot_pca_prop_var_", batch, ".pdf")),
                 exprData1 = exprData1,
                 exprData2 = exprData2)

    #pearson correlation scatter plot
    correlationPlot(exprData1 = exprData1,
                    exprData2 = exprData2,
                    batch = batch,
                    NameString = NameString)

    #boxplots before and after batch correction
    boxplot_data (expr = expr1,
                  when = "before",
                  NameString = NameString,
                  batch = batch)
    boxplot_data (expr = as.data.frame(expr2),
                  when = "after",
                  NameString = NameString,
                  batch = batch)

  }
  print("========================================================")
  print(sessionInfo())

}


#mClust turns contiguous data into discrete batches for batch correction
mclust_cluster <- function(data){
  yBIC2_Q <-mclust::mclustBIC(as.numeric(data),G=1:8)
  rs2_Q <-summary.mclustBIC(object = yBIC2_Q,data = as.numeric(data))
  tau_Q <-rs2_Q$classification
  return(tau_Q)
}
