BiocManager::install(c("affy", "GEOquery", "limma"),force =TRUE)
suppressPackageStartupMessages({
  library(affy)
  library(GEOquery)
  library(limma)
}) 

data <- getGEO(filename = "C:/GSE93303_series_matrix.txt", GSEMatrix = TRUE, getGPL = FALSE)
data = exprs(data)
data = data[, -c(9:12)] 
colnames(data) = c("DMSO_1", "DMSO_2", "DMSO_3", "DMSO_4", "RA_1", "RA_2", "RA_3", "RA_4")
head(data)

par(mfrow =  c (2,4)) # 2 rows 4 colums

hist(data[, 1],breaks = 20, main="DMSO_1",xlab ="Expression Value", ylab="Frequency")
hist(data[, 2], breaks = 20, main = "DMSO_2", xlab = "Expression value", ylab = "Frequency")
hist(data[, 3], breaks = 20, main = "DMSO_3", xlab = "Expression value", ylab = "Frequency")
hist(data[, 4], breaks = 20, main = "DMSO_4", xlab = "Expression value", ylab = "Frequency")
hist(data[, 5], breaks = 20, main = "RA_1", xlab = "Expression value", ylab = "Frequency")
hist(data[, 6], breaks = 20, main = "RA_2", xlab = "Expression value", ylab = "Frequency")
hist(data[, 7], breaks = 20, main = "RA_3", xlab = "Expression value", ylab = "Frequency")
hist(data[, 8], breaks = 20, main = "RA_4", xlab = "Expression value", ylab = "Frequency")    

par(mfrow =  c(1,1))
boxplot(data,main="Boxplots of expression values",xlab="Samples", ylab="Expression Vlues",pch =3,col = c("blue","blue","blue","blue","green","green","green","green"))

t.test(data["axo07326-r_at",c(1:4)],data["axo07326-r_at", c(5:8)])

t.test(data["axo29537-f_at", c(1:4)], data["axo29537-f_at", c(5:8)])

t.test(data["axo31319-f_s_at", c(1:4)], data["axo31319-f_s_at", c(5:8)])

t.test(data["AFFX-BioB-3_at", c(1:4)], data["AFFX-BioB-3_at", c(5:8)])

par(mfrow = c(2, 2))

barplot(data["axo07326-r_at",],main="axo07326-r_at",xlab = "Samples", ylab = "Expression value", col = c(rep("blue", 4), rep("green", 4)))

barplot(data["axo29537-f_at", ], main = "axo29537-f_at", xlab = "Samples", ylab = "Expression value", col = c(rep("blue", 4), rep("green", 4)))

barplot(data["axo31319-f_s_at", ], main = "axo31319-f_s_at", xlab = "Samples", ylab = "Expression value", col = c(rep("blue", 4), rep("green", 4)))

barplot(data["AFFX-BioB-3_at", ], main = "AFFX-BioB-3_at", xlab = "Samples", ylab = "Expression value", col = c(rep("blue", 4), rep("green", 4)))

treatment = factor(c(rep("DMSO", 4), rep("RA", 4)),levels = c("DMSO", "RA"))
design = model.matrix(~treatment)
fit = lmFit(data, design)
fit = eBayes(fit, trend = TRUE, robust = TRUE)
results = decideTests(fit)
summary(results)
topTable(fit, coef = "treatmentRA", n = 20)


resultsTable = topTable(fit, coef = "treatmentRA", n = 20080)

plot(resultsTable$logFC, -log10(resultsTable$adj.P.Val), main = "Volcano plot", xlab = "logFC", ylab = "-log10(adj.P.Val)")
abline(h = -log10(0.05), col = "red")
abline(v = c(-1, 1), col = "blue")

library(ggplot2)

ggplot(data = resultsTable, aes(x = logFC, y = -log10(adj.P.Val))) + geom_point()


resultsTable$DEG = "NO"

resultsTable$DEG[resultsTable$logFC > 0.6 & resultsTable$adj.P.Val < 0.05] <- "UP"


resultsTable$DEG[resultsTable$logFC < -0.6 & resultsTable$adj.P.Val < 0.05] <- "DOWN"

ggplot(data = resultsTable, aes(x = logFC, y = -log10(adj.P.Val), col = DEG)) + geom_point() + geom_vline(xintercept = c(-0.6, 0.6), col = "red") + geom_hline(yintercept = -log10(0.05), col = "red") + scale_color_manual(values=c("blue", "black", "red"))

plotMD(fit, coef = "treatmentRA", status = results)

sessionInfo()