#' @title Linear Regression of Principal Components
#'
#' @description Calculates Linear Regression of Principal Components 1 and 2 of the data
#' to see if the batch is correlated to the data.
#'
#' @param exprData A matrix containing gene expression data. Row names should be samples and column names should be genes.
#' @param batch.info A data frame containing the samples names and details of the batch they belong to.
#' @param batch Title of the batch the data is being corrected for.
#' @param when String indicating when the linear regression analysis is taking place
#' (before batch correction or after batch correction?). Values should be either "before" or "after".
#'
#'
#' @import stats
#'
#'
#' @return Returns the p-value associated with first ten principal components after linear regression analysis
#'
#' @export
lin_reg <- function(exprData, batch.info, batch = "Batch", NameString = "", when = "")
{
  print("===========================Linear Regression Analysis=====================")

  #calculating the principal components
  print("Calculating Principal Components...")
  pca.dat <- prcomp(exprData, center = TRUE, scale. = TRUE)

  #matching sample IDs with batches
  if (is.character(batch.info[,1])){
    id <- as.character(rownames(exprData))
    s<- match (id, batch.info[,1])

  } else if (is.numeric(batch.info[,1])){
    id <-as.numeric(rownames(exprData))
    s<- match (id, batch.info[,1])
  }

  #Adding the principal components and batch information to pca_data
  pca_data <- cbind.data.frame(batch.info[s, 2], signif(pca.dat$x[, 1:10], 5))
  colnames (pca_data)[1] <- "Batch"

  #linear regression to check which PC is associated with batch
  PCs <- colnames(pca_data[, -1])
  lm.pc <- function(x){
    print ("======================================================")
    print(paste0("Performing Linear Regression Analysis for ", x))
    pc.lm <-  lm(formula = (pca_data[,as.character(x)]) ~ as.factor(pca_data$Batch))
    print (summary(pc.lm))
    f_stat<- summary(pc.lm)$fstatistic
    p_val <- pf(f_stat[1], f_stat[2], f_stat[3], lower.tail = FALSE)
    return(p_val)
  }
  p.val <- sapply(PCs, lm.pc, simplify = TRUE)

  #plotting boxplots for first ten PCs
  plot.data <- reshape2::melt(pca_data,  id ="Batch")

  #to add p-value labels to the plot
  names (p.val) <- strsplit(names(p.val), split = ".value")
  p.val.annon <- as.character(paste0("p-value =", signif(p.val, 5)))
  p.val.dat <- vector(mode = "character", length = dim(plot.data)[1])
  index <- c((dim(plot.data)[1]/10)*c(0:9) +1)
  p.val.dat[index] <- p.val.annon

  #plotting boxplot
  date <- as.character(format(Sys.Date(), "%Y%m%d"))
  if (NameString==""){
    plotFile <- paste0(date, "_plot_", batch, "_batch_effect_boxplot_", when, "_correction.pdf")
  } else{
    plotFile <- paste0(date, "_plot_", NameString, "_", batch, "_batch_effect_boxplot_", when, "_correction.pdf")
  }
  print(paste0("Plotting boxplots showing how batches are associated with the first ten principal components to ", plotFile))
  pdf (plotFile, height = 10, width = 14)


  #boxplot to show association of batch effect with the 10 PCs
  boxplot.PC <- ggplot2::ggplot(data= plot.data, ggplot2::aes(x = Batch, y=value)) +

    ggplot2::geom_boxplot(ggplot2::aes( fill = Batch))+

    ggplot2::labs(x = "Batch", y ="Principal Components",
                  title = "Boxplots showing variation in PCs due to batch")+

    ggplot2::scale_fill_manual(values = c("#fcba03",  "#19e01c",
                                                        "#ff470f","#0fdbca",
                                                        "#ff217e","#405ce6",
                                                        "#6b6769","#b264ed"))+
    ggplot2::theme(
      #adjusting axis lines and text
      axis.line.x = ggplot2::element_line(size =0.5),
      axis.line.y = ggplot2::element_line(size =0.5),
      axis.text.x= ggplot2::element_text(size=10, colour ="black"),
      axis.text.y= ggplot2::element_text(size=10, colour ="black"),

      #Center align the title
      plot.title = ggplot2::element_text(face = "bold", hjust =0.5),

      # Adjust panel border
      panel.border = ggplot2::element_rect(fill=NA, size= 0.75),

      # Remove panel grid lines
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),

      # Remove panel background
      panel.background = ggplot2::element_blank())+
    ggplot2::geom_text(mapping = ggplot2::aes(x =2.5, y = 120),
                       label =p.val.dat,
                       size = 3) +
    ggplot2::facet_wrap(.~variable, scales = "free")

  plot(boxplot.PC)

  #plotting p-values for linear regression for Batch Variation
  p.val.data <- cbind.data.frame(as.factor(1:10), -log10(p.val))
  colnames(p.val.data)<- c("PC", "log.p.value")

  dot.plot <- ggplot2::ggplot(data = p.val.data, ggplot2::aes(x =log.p.value, y =PC)) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = -log10(0.05)), size =1,
                        linetype = 'dashed', colour = "dodgerblue3")+
    ggplot2::geom_point(size = 2.5, colour = "orange")+
    ggplot2::labs(x = "-log10(p-value)",
                  y = "Principal Components",
                  title = "Plot showing PCs associated with batch") +

    ggplot2::theme_classic()

  plot(dot.plot)
  dev.off()

  return(p.val.data)
}


