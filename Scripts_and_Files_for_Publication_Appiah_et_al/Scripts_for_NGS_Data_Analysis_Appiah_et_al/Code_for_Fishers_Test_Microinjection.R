# Run the following script before this analysis: Code_for_Microinjection_Dataset.Rmd
# Computing contribution of each condition (DMSO/EPZ)to proportion of cells in each cluster

clust <- sce_hvg.seurat$seurat_clusters
clust_pts <- data.frame(names=names(clust),cluster=clust)
clust_pts$condition <- substr(clust_pts$names, start = 0, stop = 6)
levels(clust_pts$cluster)
table_clust_pts <- table(clust_pts$condition, clust_pts$cluster)

# Perform Fishers test
datalist <- list()
data24h <- data.frame("CellsIn"=c(1,1), "CellsOut"=c(1,1))

for (j in 1:dim(table_clust_pts)[2]) {
  for (i in 1:dim(table_clust_pts)[1]) {
    data24h$CellsIn[i] <- table_clust_pts[i,j]
    data24h$CellsOut[i] <- (sum(table_clust_pts[i,]) - table_clust_pts[i,j])
    datalist[[j]] <- data24h
    rownames(datalist[[j]]) <- rownames(table_clust_pts)
  }
}

fisher_results <- list()

for (i in 1:length(datalist)) {
  fisher_results[[i]] <- fisher.test(datalist[[i]])
  
}

fisher_results
#Export results from Fishers Exact test
capture.output(fisher_results,file = "Fisher_results_aMIE14_seurat")
######
