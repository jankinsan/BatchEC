#' @title Linear Regression of Principal Components
#'
#' @description Calculates Linear Regression of Principal Components 1 and 2 of the data
#' to see if the batch is correlated to the data.
#'
#' @param exprData A matrix containing gene expression data. Row names should be samples and column names should be genes.
#' @param batch.info A data frame containing the samples names and details of the batch they belong to.
#' @param batch Title of the batch the data is being corrected for.
#'
#' @import stats
#'
#'
#' @return Returns the p-value associated with PC1 and PC2 after linear regression analysis
#'
#' @export
lin_reg <- function(exprData, batch.info, batch = "Batch")
{
  print("===========================Linear Regression Analysis=====================")

  #calculating the principal components
  print("Calculating Principal Components...")
  pca1 <- prcomp(exprData, center = TRUE, scale. = TRUE)

  #matching sample IDs with batches
  if (is.character(batch.info[,1])){
    id <- as.character(rownames(exprData))
    s<- match (id, batch.info[,1])

  } else if (is.numeric(batch.info[,1])){
    id <-as.numeric(rownames(exprData))
    s<- match (id, batch.info[,1])
  }


  #Adding the principal components and batch information to pca_data
  pca_data <- cbind.data.frame(pca1$x[, 1:2], batch.info[s, 2])
  colnames (pca_data)[3] <- "Batch"

  print("Performing Linear Regression Analysis...")

  #linear regression using PC1
  pc1.lm <- lm(formula = pca_data$PC1 ~ as.factor(pca_data$Batch))
  print(summary(pc1.lm))

  #retrieving p value
  f_stat_1<- summary(pc1.lm)$fstatistic
  p_val_PC1 <- pf(f_stat_1[1], f_stat_1[2], f_stat_1[3], lower.tail = FALSE)

  #linear regression using PC2
  pc2.lm <- lm(formula = pca_data$PC2 ~ as.factor(pca_data$Batch))
  print(summary(pc2.lm))

  #retrieving p value
  f_stat_2<- summary(pc2.lm)$fstatistic
  p_val_PC2 <- pf(f_stat_2[1], f_stat_2[2], f_stat_2[3],
                  lower.tail = FALSE)

  print(paste0("P-value for ", batch, " from PC1: ", p_val_PC1, " from PC2: ", p_val_PC2))
  p_val <- c(p_val_PC1, p_val_PC2)

  return(p_val)
}



