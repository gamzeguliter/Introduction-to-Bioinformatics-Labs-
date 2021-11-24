#quýz 1 ELÝF GAMZE GÜLÝTER
#QUESTION 1
quiz_data <- read.csv("C:/Users/gamze.guliter-ug/Downloads/quiz_exp_file.csv", sep=",", header=TRUE)
quiz_data
class(quiz_data)
colnames(quiz_data)
rownames(quiz_data)
length(quiz_data)
dim(quiz_data)
normal = quiz_data[,2:11]
normal
tumor = quiz_data[,12:21]
tumor[4,] #Mek2
normal[9,] #MTORC2

#QUESTION 2
#NOTE : the mean() function gave errors therefor ý calculated mean myself 

mean_n = normal[9,1] + normal[9,2]+ normal[9,3]+normal[9,4]+normal[9,5]+normal[9,6]+normal[9,7]+normal[9,8] + normal[9,9] + normal[9,10]
mean_n = mean_n / 10
mean_n

help(mean)
max_n = max(normal[9,])
max_n
median_n = median(normal[9,])
min_n = min(normal[9,])
stand_n = sd(normal[9,])

mean_t = tumor[4,1] + tumor[4,2]+ tumor[4,3]+tumor[4,4]+tumor[4,5]+tumor[4,6]+tumor[4,7]+tumor[4,8] + tumor[4,9] + tumor[4,10]
mean_t= mean_t / 10
mean_t
max_t = max(tumor[4,])
max_t
median_t = median(tumor[4,])
min_t = min(tumor[4,])
stand_t = sd(tumor[4,])

info_t= c(mean_t, max_t, median_t, min_t, stand_t)
info_t


info_n= c(mean_n, max_n, median_n, min_n, stand_n)
info_n

normal <-  rbind(normal, info_n) 
normal


tumor <-  rbind(tumor, info_t) 
tumor
#QUESTION 3
mean_difference = mean_n - mean_t
mean_difference

#QUESTION 4
#gene expression is higher for the PI3K pathaway components. the pathaway starting from the AKT to mTROC1 has the hihgest variablity