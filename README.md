# BatchEC
An R package that evaluates batch effects in gene expression data and corrects batches using the ComBat algorithm (available in sva package).

## Usage
install_github("jankinsan/BatchEC")


library(BatchEC)

***#if input is to be given as tab delimited text files***
batch_eval_cor(
  exprFile = "~filepath.txt",
  batchFile = "~filepath2.txt",
  batch,
  NameString = "",
  discrete.batch = TRUE)

***#if input objects already exist in the R environment***
beacon(expr1, 
batch.info, 
batch,
NameString = "",
discrete.batch = TRUE)
	where expr1 is gene expression data; rows should be genes and columns should be samples.

## Functions
**batch_eval_cor**	Batch effects evaluation and correction from files. Both the files (containing the Expression data and the batch information should be tab-delimited text files.)


**beacon**	        	Batch effects evaluation and correction from objects loaded in the workspace.


**boxplot_data**		Plots boxplots showing Gene Expression across samples


**ComBat_data**		Batch Correction using ComBat


**correlationPlot**   	Correlation Scatter Plot


**kmeans_PCA**      	PCA plot with k-means information


**lin_reg**            	Linear Regression of Principal Components


**pca_batch**           	PCA Plot with batch information


**pca_prop_var**        	PCA Proportion of Variation
