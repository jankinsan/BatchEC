#' @title Batch Correction using ComBat
#'
#' @description Corrects the gene expression data for a particular batch using
#' ComBat algorithm
#'
#'
#' @param expr gene expression data to be corrected for batch effects. Rows
#'  should be genes and columns should be samples.
#' @param batch.info A data frame containing batch information corresponding to
#' the gene expression dataset,the first column should be sample names and the
#' second column should contain the batch the sample belongs to.
#' @param NameString  String that will be added to the names of all output files. Default=NA.
#' @param batch The title of the batch for which you want to evaluate and do the
#'  correction. Default = "Batch"
#'
#' @import utils
#'
#'
#' @return Returns a matrix containing the batch corrected values for gene expression.
#'
#' @export
ComBat_data <- function(expr, batch.info, batch = "Batch", NameString = "")
{
  print ("===========================Batch Effects Adjustment using ComBat=====================")
  #matching IDs for ComBat
  if (is.character(batch.info[,1])){
    match.id <- match (as.character(colnames(expr)), batch.info[,1])
    batch.id <- batch.info[match.id, 2]

  } else if (is.numeric(batch.info[,1])){
    match.id <- match (as.numeric(colnames(expr)), batch.info[,1])
    batch.id <- batch.info[match.id, 2]
  }



  print("Performing batch correction using ComBat...")
  batch_corrected <- sva::ComBat(dat=expr,
                                 batch = batch.id,
                                 mod=NULL,
                                 par.prior=TRUE,
                                 prior.plots=TRUE)

  date <- as.character(format(Sys.Date(), "%Y%m%d"))
  if(NameString==""){
    outfile  <- outFile <- paste0(date, "_data_", "batch_corrected_", batch, ".txt")
  } else {
    outFile <- paste0(date, "_data_", NameString, "_batch_corrected_", batch, ".txt")
  }

  write.table(batch_corrected, outFile, quote=FALSE, sep="\t")
  print(paste("Batch corrected data written to file... ", outFile, sep=""))
  return (batch_corrected)
}
