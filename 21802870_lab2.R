genes = c("PTEN", "TP53", "HER2", "GR", "PGR", "MDM4")
expression =c("3.89", "-5.90", "18", "38.48", "-0.52", "22.5")
genes
expression
expression_data  = cbind(genes,expression)
expression_data
colnames(expression_data) 


colnames(expression_data) = c( "gene_ids" , "log2CPM")
expression_data_compare = expression_data[,"log2CPM"] >= 2 
expression_data_compare


expression_data_filtered = expression_data[expression_data_compare,]
expression_data_filtered 

class(expression_data_filtered )
orthology_data = read.csv("D:/Users/gamze.guliter-ug/Downloads/worksheet_sample_data.csv", sep=",", header=TRUE)

orthology_data
comp = orthology_data[,3] > 0 
comp 
orthology_data_filtered =orthology_data[comp,]
orthology_data_filtered 


write.csv(orthology_data_filtered, "orthology_data_filtered.xlsx")