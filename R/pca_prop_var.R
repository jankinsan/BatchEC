#' @title PCA Proportion of Variation
#'
#' @description Calculates and plots the proportion of variation for twenty-five Principal Components of
#' expression data before and after batch correction
#'
#' @param batch.title The tile of the batch
#' @param plotFile The PCA proportion of variation will be saved to this file
#' @param exprData1 A matrix containing the gene Expression data before batch correction
#' @param exprData2 A matrix containing the gene Expression data after batch correction
#'
#' @import grDevices
#' @import stats
#'
#' @export
pca_prop_var <- function(batch.title, plotFile, exprData1, exprData2){

  print("===========================Plotting PCA Proportion of Variation=====================")

  #calculating the principal components
  print("Calculating Principal Components...")
  pca1 <- prcomp(exprData1, center = TRUE, scale. = TRUE)
  pca2 <- prcomp(exprData2, center = TRUE, scale. = TRUE)

  before_correction <- summary(pca1)$importance[2,]
  after_correction <- summary(pca2)$importance[2,]

  #making a dataframe from proportion of variation values
  s <- data.frame(as.factor(names(before_correction)),
                  before_correction,
                  after_correction,
                  row.names = 1:dim(exprData1)[1])
  colnames(s)[1] <- "PC"
  s$PC <- factor(s$PC, levels = s$PC[order(s$before_correction, decreasing=TRUE)])

  #plotting Proportion of Variation
  print("Plotting Proportion of Variation")
  date <- as.numeric(format(Sys.Date(), "%Y%m%d"))
  pdf(paste0(date, "_", plotFile))
  prop_var_plot <- ggplot2::ggplot(data=s[1:15,])+

    ggplot2::geom_point(show.legend = TRUE, size =3,
                        ggplot2::aes(PC, before_correction, colour ="dodgerblue3", group = 1))+

    ggplot2::geom_path(size=1.5, ggplot2::aes(PC, before_correction, group=1, colour ="dodgerblue3"))+

    ggplot2::geom_point(show.legend = TRUE, size =3,
                        ggplot2::aes(PC, after_correction, colour ="orange", group = 2))+

    ggplot2::geom_path(size=1.5, ggplot2::aes(PC, after_correction, group=2, colour ="orange"))+


    ggplot2::scale_colour_manual(values =c("dodgerblue3" ,"orange"),
                                 labels = c('Before Correction','After Correction'))+

    ggplot2::labs(x="Principal Components", y="Proportion of Variation",
                  title = paste0("Proportion of Variation for ",batch.title))+

    ggplot2::theme(

      #adjusting axis lines and text
      axis.title = ggplot2::element_text(size = 15),
      axis.line.x = ggplot2::element_line(size =0.75),
      axis.line.y = ggplot2::element_line(size =0.75),
      axis.text.x = ggplot2::element_text(angle = 90, size=15, colour ="black"),
      axis.text.y = ggplot2::element_text(size=15, colour ="black"),

      #Center align the title
      plot.title = ggplot2::element_text(face = "bold", hjust =0.5, size =20),

      #modify legend
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 15),
      legend.position = c(0.8, 0.8),

      # Adjust panel border
      panel.border = ggplot2::element_blank(),

      # Remove panel grid lines
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),

      # Remove panel background
      panel.background = ggplot2::element_blank()

    )


  plot(prop_var_plot)
  dev.off()

}
