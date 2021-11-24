BiocManager::install(c("affy", "GEOquery", "edgeR"), force=TRUE)
suppressPackageStartupMessages({
  library(affy)
  library(GEOquery)
  library(edgeR)
}) 

axolotl_data <- getGEO(filename = "C:/Users/gamze.guliter-ug/Downloads/GSE93303_series_matrix.txt/GSE93303_series_matrix.txt", GSEMatrix = TRUE, getGPL = FALSE)
axolotl_data <- exprs(axolotl_data)
head(axolotl_data)

colnames(axolotl_data) <- c("DMSO_1", "DMSO_2", "DMSO_3", "DMSO_4",
                            "RA_1", "RA_2", "RA_3", "RA_4",
                            "LE135_1","LE135_2","LE135_3","LE135_4")
head(axolotl_data)
dim(axolotl_data) 
GSE93303_gene_names <- read.table("C:/Users/gamze.guliter-ug/Desktop/GSE93303_gene_names.csv", sep=";", header=TRUE)
head(GSE93303_gene_names) 
dim(GSE93303_gene_names)  
GSE93303_gene_names <- GSE93303_gene_names[!apply(GSE93303_gene_names == "" , 1, any), ] 
head(GSE93303_gene_names)
dim(GSE93303_gene_names)

GSE93303_gene_names <- GSE93303_gene_names[!duplicated(GSE93303_gene_names$GeneSymbol),]
dim(GSE93303_gene_names)

axolotl_data <- merge(GSE93303_gene_names, axolotl_data, by.x="ID", by.y=0) 
head(axolotl_data)

rownames(axolotl_data) <- axolotl_data$GeneSymbol
axolotl_data <- axolotl_data[,-c(1:2)]
set.seed(5432) 
hist(distribution, main="Normal distribution, (n=100, mean=12, sd=1)", xlab="values", col="orange") 

logCPM <- cpm(axolotl_data, log=TRUE)
variance <- apply(axolotl_data, 1, var)
head(variance)

variance_selected <- names(sort(variance, decreasing=TRUE))[1:500]
head(variance_selected)

axolotl_data_top_500 <- axolotl_data[variance_selected,]
dim(axolotl_data_top_500)
head(axolotl_data_top_500)

compute_variance <- function(x) sd(x)^2
variance2 <- apply(axolotl_data, 1, compute_variance)
head(variance2)

compute_cv <- function(x) sd(x) / mean(x)
cv <- apply(axolotl_data, 1, compute_cv)
head(cv)