# BatchEC
An R package that evaluates batch effects in gene expression data and corrects batches using the ComBat algorithm (available in sva package).

# Usage
install_github("jankinsan/BatchEC")


library(BatchEC)

# Functions
batch_eval_cor   Batch effects evaluation and correction


boxplot_data      Plots boxplots showing Gene Expression across samples


ComBat_data       Batch Correction using ComBat


correlationPlot   Correlation Scatter Plot


kmeans_PCA        PCA plot with k-means information


lin_reg            Linear Regression of Principal Components


pca_batch           PCA Plot with batch information


pca_prop_var        PCA Proportion of Variation
