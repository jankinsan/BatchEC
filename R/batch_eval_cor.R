#' @title Batch Evaluation and Correction
#'
#'
#' @description Evaluates if the batch effects with a gene expression dataset using linear
#' regress and if the batch is associated with the data, the data is batch-corrected
#' using the ComBat Algorithm.
#'
#' Please set your working directory before you call the function.
#' The directory should contain the input files.
#'
#'
#' @param exprFile   File name of tab-delimited .txt file that contains the expression data. In the file, rows
#' should be genes and columns should be samples.
#' @param batchFile Path to the tab-delimited .txt file that contains the batch information.
#' The first column should be sample names.
#' @param batch The title of the batch for which you want to evaluate and do the
#'  correction. Default = "Batch"
#' @param NameString  String that will be added in all output filenames. Default=NA.
#' @param discrete.batch Logical value indicating whether the samples are already
#' grouped in discrete batches. If the value is FALSE, contiguous batch information is
#' clustered into discrete batches. Useful for clustering batch variables that are
#' contiguous, for example, used reads and useful reads in mClust. Default = TRUE
#'
#' @import utils
#' @import mclust
#'
#' @examples
#' \dontrun{
#' batch_eval_cor(exprFile = "~/exprData.txt",
#'                batchFile = "~/batchData.txt",
#'                Batch = "Batch",
#'                discrete.batch = TRUE)
#' }
#'
#'
#' @export
batch_eval_cor <- function(exprFile, batchFile, batch, NameString = "", discrete.batch = TRUE){

  #reading expression data from file
  print (paste0("Reading gene expression data from ", exprFile))
  expr1 <- read.table(exprFile,
                    header = TRUE,
                    stringsAsFactors = FALSE,
                    row.names = 1,
                    check.names=FALSE)

  #reading batch information
  print (paste0("Reading batch data from ", batchFile))
  batch.info <- read.table(batchFile,
                           header = TRUE,
                           stringsAsFactors = FALSE,
                           sep = "\t")
  #matching batch
  mat <- match(batch, colnames(batch.info))
  batch.info <- batch.info [,c(1, mat)]

  if (discrete.batch== FALSE){
    batch.info[,2] <- mclust_cluster(batch.info[,2])

    #writing the batch cluster information to file
    date <- as.character(format(Sys.Date(), "%Y%m%d"))
    clusterInfoFile <- paste0(date, "_", NameString, "_", batch, "_batch_info_mClust.txt")
    print(paste0("Writing the batch information from mClust for ", batch, "to: ", clusterInfoFile))
    write.table(batch.info,
                file = clusterInfoFile,
                sep = "\t",
                col.names = TRUE,
                row.names = FALSE)
  }

  exprData1 <- t(expr1)

  #linear regression before batch correction
  p_val_before <- lin_reg(exprData=exprData1,
                          batch.info = batch.info,
                          batch=batch)

  #checking if any of the p values are less than 0.05
  if(p_val_before[1] >0.05 && p_val_before[2] >0.05){
    stop("Halting execution since the batch specified is not associated with data...")

  }else{

    print("Batch is associated with the data...")

    #pca with batch before correction
    pca_batch (exprData = exprData1,
               batch.info= batch.info,
               batch= batch,
               NameString = NameString,
               when = "before_ComBat_correction")

    #pca with batch and kmeans before batch correction
    k_before <- kmeans_PCA(exprData = exprData1,
                          batch.info = batch.info,
                          batch = batch,
                          NameString = NameString,
                          when = "before_correction")

    #Batch Correction using ComBat
    expr2 <- ComBat_data (exprData = expr1,
                         batch.info= batch.info,
                         batch = batch,
                         NameString = NameString)
    exprData2 <- t(expr2)

    #linear regression after batch correction
    p_val_after <- lin_reg(exprData=exprData2, batch.info = batch.info, batch=batch)

    #pca with batch after correction
    pca_batch (exprData = exprData2,
               batch.info= batch.info,
               batch= batch,
               NameString = NameString,
               when = "after_ComBat_correction")

    #pca with batch and k-means after correction
    k_after <- kmeans_PCA(exprData = exprData2,
                          batch.info = batch.info,
                          batch = batch,
                          NameString = NameString,
                          when = "after_correction")

    #ensuring biological variability is not removed by checking k before and after correction
    if (k_before != k_after){
      stop("Optimal number of clusters differ before and after correction, biological variability may have been removed by ComBat")
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
                    fileName = NameString)

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

}


#mClust turns contiguous data into discrete batches for batch correction
mclust_cluster <- function(data){
  yBIC2_Q <-mclust::mclustBIC(as.numeric(data),G=1:8)
  rs2_Q <-summary.mclustBIC(object = yBIC2_Q,data = as.numeric(data))
  tau_Q <-rs2_Q$classification
  return(tau_Q)
}
