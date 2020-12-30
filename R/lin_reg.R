#' @title Linear Regression of Principal Components
#'
#' @description Calculates Linear Regression of Principal Components 1 and 2 of the data
#' to see if the batch is correlated to the data.
#'
#'
#' @param data A matrix containing gene expression data. Row names should be
#' samples and column names should be genes.
#' @param batch.info A data frame containing the samples names and details of the
#' batch they belong to.
#' @param batch Title of the batch the data is being corrected for.
#'
#' @import stats
#'
#'
#' @return Returns the p-value after linear regression analysis if the batch is
#' associated with data.
#'
#' @export
lin_reg <- function(data, batch.info, batch = "Batch")
{
  print("===========================Linear Regression Analysis=====================")
  #calculating the principal components
  print("Calculating Principal Components...")
  pca1 <- prcomp(data, center = TRUE, scale. = TRUE)

  #matching sample IDs with batches
  id <- as.numeric(rownames(data))
  s<- match (id, batch.info[,1])

  #Adding the principal components and batch information to pca_data
  pca_data <- cbind.data.frame(pca1$x[, 1:2], batch.info[s, 2])
  colnames (pca_data)[3] <- "Batch"

  print("Performing Linear Regression Analysis...")


  #plots to see the gene expression as a function of batch
  plot(as.factor(pca_data$Batch), pca_data$PC1)

  #linear regression using PC1
  pc.lm1 <- lm(formula = pca_data$PC1 ~ as.factor(pca_data$Batch))
  print(summary(pc.lm1))

  #retrieving p value
  f_stat_1<- summary(pc.lm1)$fstatistic
  p_val_PC1 <- pf(f_stat_1[1], f_stat_1[2], f_stat_1[3], lower.tail = FALSE)

  #linear regression using PC2
  pc.lm2 <- lm(formula = pca_data$PC2 ~ as.factor(pca_data$Batch) )
  print(summary(pc.lm2))

  #retrieving p value
  f_stat_2<- summary(pc.lm2)$fstatistic
  p_val_PC2 <- pf(f_stat_2[1], f_stat_2[2], f_stat_2[3],
                  lower.tail = FALSE)

  print(paste0("P-value for ", batch, " from PC1: ", p_val_PC1, " from PC2: ", p_val_PC2))
  p_val <- c(p_val_PC1, p_val_PC2)

  #checking if both the p values are less than 0.05
  max_p_val<- which.max(p_val)
  if(max_p_val <0.05){
    print("Batch is associated with the data")
    return(p_val)
  } else{
    print("Batch is not associated with the data")
    return(NA)
    }
}



