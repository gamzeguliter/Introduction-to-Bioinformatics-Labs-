BiocManager::install(c("affy","GEOquery"))
suppressPackageStartupMessages({
  library(affy)
  library(GEOquery)
}) 
for (i in 1:5) {
  print(i)
}
weekdays = c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday")

for (day in weekdays) {
  print(day)
}

#C:\Users\gamze.guliter-ug\Downloads\GSE93303_series_matrix (1).txtdata <- getGEO(filename = "GSE93303_series_matrix.txt", GSEMatrix = TRUE, getGPL = FALSE)
data <- getGEO(filename = "C:/Users/gamze.guliter-ug/Desktop/GSE93303_series_matrix.txt/GSE93303_series_matrix.txt", GSEMatrix = TRUE, getGPL = FALSE)
data = exprs(data)
data = data[,-c(9:12)]
head(data)
j = 0
 while(j <=7){
   print(j)
   j = j +1
 }
mat = cbind(c(3,0,3,3),c(3,0,0,0),c(3,0,0,3),c(1,1,0,0),c(1,1,1,0),c(1,1,1,0))
mat
result_1 =apply(X=mat,MARGIN =1 , FUN = sum)
result_1


result_2 =apply(X=mat,MARGIN =2 , FUN = sum)
result_2

sd_values = apply(X = data, MARGIN = 1, FUN = sd)
head(sd_values)


data.df = as.data.frame(data)
data.df$mean_description = ""
head(data.df)

for (i in 1:dim(data.df)[1]) {
  if (mean(as.numeric(data.df[i, 1:8])) <= 4) {
    data.df$mean_description[i] <- "LOW"
  } else if (mean(as.numeric(data.df[i, 1:8])) > 4 & mean(as.numeric(data.df[i, 1:8])) < 10) {
    data.df$mean_description[i] <- "AVERAGE"
  } else {
    data.df$mean_description[i] <- "HIGH"
  }
}

head(data.df)


subset.df = data.frame()

index = 1

while (dim(subset.df)[1] < 1000) {
  if (mean(as.numeric(data.df[index, 1:4])) > mean(as.numeric(data.df[index, 5:8]))) {
    subset.df <- rbind.data.frame(subset.df, data.df[index, ])
  }
  
  index = index + 1
}

head(subset.df)
sessionInfo()