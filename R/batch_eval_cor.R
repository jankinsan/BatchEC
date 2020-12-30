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
#' @param Batch The title of the batch for which you want to evaluate and do the
#'  correction. Default = "Batch"
#' @param batchFile Path to the tab-delimited .txt file that contains the batch information.
#' The first column should be sample names.
#' @param NameString  String that will be added in all output filenames. Default=NA.
#' @param discrete.batch Logical value indicating whether the samples are already
#' grouped in discrete batches. If the value is FALSE, contiguous batch information is
#' clustered into discrete batches. Useful for clustering batch variables that are
#' contiguous, for example, used reads and useful reads in mClust. Default = TRUE
#'
#' @import utils
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
batch_eval_cor <- function(exprFile, batchFile, Batch, NameString = "", discrete.batch = TRUE){

  #reading expression data from files
  print (paste0("Reading gene expression data from ", exprFile))
  exp1 <- read.table(exprFile,
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
  mat <- match(Batch, colnames(batch.info))
  batch.info <- batch.info [,c(1, mat)]

  if (discrete.batch== FALSE){
    batch.info[,2] <- mclust_cluster(batch.info[,2])
  }

  data1 <- t(exp1)

  #linear regression
  p_val_before <- lin_reg(data=data1, batch.info = batch.info, batch=Batch)

  if(is.na(p_val_before)){
    print("Halting execution since the batch specified is not associated with data...")

  } else {

    #pca with batch before correction
    pca_batch (data = data1,
               batch.info= batch.info,
               batch= Batch,
               plotFile = ifelse (NameString =="",
                                  paste0("plot_", Batch, "_pca_before_batch_correction.pdf"),
                                  paste0("plot_", NameString, "_", Batch, "_pca_before_batch_correction.pdf")))

    #pca with batch and kmeans before batch correction
    kmeans_PCA(exprData = data1,
               batch.info = batch.info,
               batch = Batch,
               NameString = NameString,
               when = "before_correction")

    #Batch Correction using ComBat
    exp2 <- ComBat_data (exprData = exp1,
                         batch.info= batch.info,
                         batch = Batch,
                         NameString = NameString)
    data2 <- t(exp2)

    #pca with batch after correction
    pca_batch (data = data2,
               batch.info= batch.info,
               batch= Batch,
               plotFile = ifelse( NameString =="",
                                  paste0("plot_", Batch, "_pca_after_batch_correction.pdf"),
                                  paste0("plot_", NameString,"_", Batch, "_pca_after_batch_correction.pdf")))

    #pca with batch and k-means after correction
    kmeans_PCA(exprData = data2,
               batch.info = batch.info,
               batch = Batch,
               NameString = NameString,
               when = "after_correction")

    #PCA Proportion of Variation
    pca_prop_var(batch.title = Batch,
                 plotFile = ifelse( NameString == "",
                                    paste0("plot_pca_prop_var_", Batch, ".pdf"),
                                    paste0(NameString, "_plot_pca_prop_var_", Batch, ".pdf")),
                 data1 = data1,
                 data2 = data2)

    #pearson correlation
    correlationPlot(data1 = data1,
             data2 = data2,
             batch = Batch,
             fileName = NameString)

    #boxplots before and after batch correction
    boxplot_data (expr1 = exp1, expr2 = exp2, NameString = NameString,
                  batch = Batch)
  }

  print(sessionInfo())
}


#mClust turns contiguous data into discrete batches for batch correction
mclust_cluster <- function(data){
  yBIC2_Q <-mclust::mclustBIC(as.numeric(data),G=1:8)
  rs2_Q <-summary(yBIC2_Q,as.numeric(data))
  tau_Q <-rs2_Q$classification
  return(tau_Q)
}
