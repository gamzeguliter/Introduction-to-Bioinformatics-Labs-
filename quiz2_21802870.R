BiocManager::install(c("affy", "GEOquery", "limma"))
suppressPackageStartupMessages({
  library(affy)
  library(GEOquery)
  library(limma)
}) 

#Question 1
data <- read.csv("C:/Users/gamze.guliter-ug/Downloads/count_data.csv", sep=",", header=TRUE)
data

length(data)

par(mfrow = c(2, 3)) # sets up a graphical container with 2 rows and 4 columns

hist(data[, 2], breaks = 20, main = "Control1", xlab = "Expression value", ylab = "Frequency")

hist(data[, 3], breaks = 20, main = "Control2", xlab = "Expression value", ylab = "Frequency")

hist(data[, 4], breaks = 20, main = "Control3", xlab = "Expression value", ylab = "Frequency")

hist(data[, 5], breaks = 20, main = "E21", xlab = "Expression value", ylab = "Frequency")

hist(data[, 6], breaks = 20, main = "E22", xlab = "Expression value", ylab = "Frequency")

hist(data[, 7], breaks = 20, main = "E23", xlab = "Expression value", ylab = "Frequency")

par(mfrow = c(1, 1))

boxplot(data[,2:7], main = "Boxplots of expression values", xlab = "Samples", ylab = "Expression values", pch = 3, col = c("blue","blue","blue","green","green","green"))
#Question 2
mean_info = data[,1:3]
head(mean_info)
colnames(mean_info)=c("gene_names","mean_Control","mean_E2")
head(mean_info)
#means of control

for (j in 2:dim(data)[1] ){
  for (i in 2:4 ){
  row_average <- mean(data[j,2:4])
  print(paste("col", i, "mean:", row_average))
}
}
#means of other

for (j in 2:dim(data)[1] ){
  for (i in 2:4 ){
    row_average <- mean(data[j,5:7])
    print(paste("col", i, "mean:", row_average))
  }
}
# Question 4
function(df,genes="all",method ="limma",control_label,treatment_label){
  if(method=="limma"){
    #perform limma version
  }  
  if(method=="ttest"){
    #perform ttest version
    if(genes==all){
      for(i in 1:m)
        t.test(df[i,], data["AFFX-BioB-3_at", c(5:8)])
    }
     
    else{    
      for(i in 1:m)
        t.test(df[,genes], data["AFFX-BioB-3_at", c(5:8)])
      
    }
    }
}