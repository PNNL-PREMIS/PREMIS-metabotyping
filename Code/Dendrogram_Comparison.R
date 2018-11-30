library(dendextend)
#Compare Albert's above and belowground metabolite clustering results
setwd("~/Documents/iPASS/Manuscripts/GoBrachy MS1/Code/")
above <- read.csv("DistanceMatrix_Aboveground_GenotypeAverages.csv")
rownames(above) <- colnames(above)[-1]
names <- rownames(above)
above <- above[,-1]
above <- data.matrix(above)
above <- above[lower.tri(above)]
attr(above,"Size") <- 30
attr(above,"Labels") <- names


below <- read.csv("DistanceMatrix_Belowground_GenotypeAverages.csv")
rownames(below) <- colnames(below)[-1]
names <- rownames(below)
below <- below[,-1]
below <- data.matrix(below)
below <- below[lower.tri(below)]
attr(below,"Size") <- 30
attr(below,"Labels") <- names


#Dendograms
above_clust <- hclust(above,method = "ward.D")
plot(above_clust)
below_clust <- hclust(below,method = "ward.D")
plot(below_clust)

#Compare dendograms with Baker's gamma
obs_gamma <- cor_bakers_gamma(above_clust,below_clust)

#Assess significance using permutation test
B <- 1000
gamma_dist <- rep(NA,B)
above_clust_perm <- above_clust
below_clust_perm <- below_clust

for(i in 1:B){
  above_clust_perm$labels <- sample(above_clust_perm$labels)
  below_clust_perm$labels <- sample(below_clust_perm$labels)
  gamma_dist[i] <- cor_bakers_gamma(above_clust_perm,below_clust_perm)
  
}

hist(gamma_dist,xlim=c(-1,1));abline(v=obs_gamma)
