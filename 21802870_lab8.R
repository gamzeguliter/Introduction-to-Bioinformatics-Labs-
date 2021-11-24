BiocManager::install("ape")
install.packages("phangorn")

library("ape")
library("phangorn")

sim_matrix <- matrix(c(1.0, 0.5, 0.3,
                       0.5, 1.0, 0.4,
                       0.3, 0.4, 1.0),
                     nrow = 3)

row.names(sim_matrix) <- c("G", "T", "M")
colnames(sim_matrix) <- c("G", "T", "M")


dissim_matrix <- 1 - sim_matrix

dist_matrix <- as.dist(dissim_matrix)

nj_tree <- ape::nj(dist_matrix)
plot(nj_tree)
plot(nj_tree, type = "unrooted")
upgma_tree <- phangorn::upgma(dist_matrix)
plot(upgma_tree)

plot(nj_tree, type = "unrooted")
upgma_tree <- phangorn::upgma(dist_matrix)
plot(upgma_tree)
par(mfrow = c(1, 2))
plot(nj_tree)
plot(upgma_tree)five_sim_matrix <- matrix(c(1.0, 0.0, 0.0, 0.0, 0.0,        
                                            0.9, 1.0, 0.0, 0.0, 0.0,                
                                            0.8, 0.7, 1.0, 0.0, 0.0,      
                                            0.5, 0.4, 0.3, 1.0, 0.0,        
                                            0.3, 0.2, 0.1, 0.8, 1.0),
                                          nrow = 5, byrow = T)

row.names(five_sim_matrix) <- c("ME", "B", "G", "T", "MW")
colnames(five_sim_matrix) <- c("ME", "B", "G", "T", "MW")
five_dissim_matrix <- 1 - five_sim_matrix

five_dist_matrix <- as.dist(five_dissim_matrix)

five_nj_tree <- ape::nj(five_dist_matrix)
five_upgma_tree <- phangorn::upgma(five_dist_matrix)

par(mfrow = c(1, 2))
plot(five_nj_tree)
plot(five_nj_tree, "unrooted")

plot(five_upgma_tree)
plot(five_upgma_tree, "unrooted")


par(mfrow = c(1, 2))
plot(five_nj_tree)
plot(five_upgma_tree)



#Q2) Compare the UPGMA and NJ trees with respect to their a) topology and b) branch lengths and c) root and outgroup.
# first of all nj tree has more nodes and branches and has a greater height. As branch length NJ tree again, has longer branch lengths. 
# I think both of them has root. As outgroup in NJ tree, outgorup is B. In ourgroup of the upgma tree, it is B or ME . 
