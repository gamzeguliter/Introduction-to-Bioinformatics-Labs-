
BiocManager::install(c("affy", "GEOquery", "limma"),force =TRUE)
suppressPackageStartupMessages({
  library(affy)
  library(GEOquery)
  library(limma)
}) 

data = getGEO(filename = "C:/GSE79299_series_matrix.txt", GSEMatrix = TRUE, getGPL = FALSE)
data = exprs(data)

data = data[, -c(7:13)] 
head(data)

colnames(data) = c("uninjured_tissue_1", "uninjured_tissue_2", "uninjured_tissue_3", "2_day_post_injury_1", "2_day_post_injury_2", "2_day_post_injury_3")

par(mfrow =  c (2,3)) # 2 rows 3 columns
hist(data[, 1],breaks = 20, main="uninjured_tissue_1",xlab ="Expression Value", ylab="Frequency")
hist(data[, 2], breaks = 20, main = "uninjured_tissue_2", xlab = "Expression value", ylab = "Frequency")
hist(data[, 3], breaks = 20, main = "uninjured_tissue_3", xlab = "Expression value", ylab = "Frequency")
hist(data[, 4], breaks = 20, main = "2_day_post_injury_1", xlab = "Expression value", ylab = "Frequency")
hist(data[, 5], breaks = 20, main = "2_day_post_injury_2", xlab = "Expression value", ylab = "Frequency")
hist(data[, 6], breaks = 20, main = "2_day_post_injury_3", xlab = "Expression value", ylab = "Frequency")

#it looks like normal distribution but the expression values between 3-5 looks outliner

par(mfrow =  c(1,1))
boxplot(data,main="Boxplots of expression values",xlab="Samples", ylab="Expression Vlues",pch =3,col = c("blue","blue","blue","blue","green","green","green","green"))

t.test(data["axo09654-f_at",c(1:3)],data["axo09654-f_at", c(4:6)])

t.test(data["axo11065-f_at",c(1:3)],data["axo11065-f_at", c(4:6)])

t.test(data["axo17266-f_at",c(1:3)],data["axo17266-f_at", c(4:6)])

t.test(data["axo15862-f_at",c(1:3)],data["axo15862-f_at", c(4:6)])

t.test(data["axo20884-f_s_at",c(1:3)],data["axo20884-f_s_at", c(4:6)])

t.test(data["axo06507-r_at",c(1:3)],data["axo06507-r_at", c(4:6)])

par(mfrow =  c(2,3))

barplot(data["axo09654-f_at",],main="axo09654-f_at",xlab = "Samples", ylab = "Expression value", col = c(rep("blue", 3), rep("green", 3)))

barplot(data["axo11065-f_at", ], main = "axo11065-f_at", xlab = "Samples", ylab = "Expression value", col = c(rep("blue", 3), rep("green", 3)))

barplot(data["axo17266-f_at", ], main = "axo17266-f_at", xlab = "Samples", ylab = "Expression value", col = c(rep("blue", 3), rep("green", 3)))

barplot(data["axo15862-f_at", ], main = "axo15862-f_at", xlab = "Samples", ylab = "Expression value", col = c(rep("blue", 3), rep("green", 3)))

barplot(data["axo20884-f_s_at", ], main = "axo20884-f_s_at", xlab = "Samples", ylab = "Expression value", col = c(rep("blue", 3), rep("green", 3)))

barplot(data["axo06507-r_at", ], main = "axo06507-r_at", xlab = "Samples", ylab = "Expression value", col = c(rep("blue", 3), rep("green", 3)))

#Analyze the dataset using the limma package. List the top 10 most significantly differentially expressed genes. Write out the results in a .csv file and upload it alongside your R script.