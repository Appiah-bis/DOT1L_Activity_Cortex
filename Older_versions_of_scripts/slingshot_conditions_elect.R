# Using Slingshot for pseudotime analysis 
#Slingshot asnalysis for DMSO and EPZ cells separately

library(slingshot)

sce_hvg.seurat_DMSO_2 <- subset(sce_hvg.seurat_DMSO, idents = c("0", "1", "2", "3", "4"))

sce_hvg.slingshot_DMSO <- as.SingleCellExperiment(sce_hvg.seurat_DMSO_2)



## Gaussian mixture model for clustering  


#cl1 <- Mclust(reducedDim(sce_hvg.slingshot)[,1:2])$classification
cl1_dmso <- Idents(sce_hvg.seurat_DMSO_2)
colData(sce_hvg.slingshot_DMSO)$GMM <- cl1_dmso

plot(reducedDim(sce_hvg.slingshot_DMSO)[,1:2], col = brewer.pal(6,"Dark2")[cl1_dmso], pch=16, asp = 1)



sce_hvg.slingshot_DMSO <- slingshot::slingshot(sce_hvg.slingshot_DMSO, clusterLabels = 'GMM', reducedDim = 'UMAP')


#saveRDS(file="sce_hvg.slingshot", sce_hvg.slingshot)
#sce_hvg.slingshot <- readRDS(file="sce_hvg.slingshot")

#summary(sce_hvg.slingshot_DMSO$slingPseudotime_1)


colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce_hvg.slingshot_DMSO$slingPseudotime_1, breaks=100)]

plot(reducedDims(sce_hvg.slingshot_DMSO)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce_hvg.slingshot_DMSO), lwd=2, col='black')



plot(reducedDims(sce_hvg.slingshot_DMSO)$PCA, col = brewer.pal(9,'Set1')[sce_hvg.slingshot_DMSO$seurat_clusters], pch=16, asp = 1)
lines(SlingshotDataSet(sce_hvg.slingshot_DMSO), lwd=2, type = 'lineages', col = 'black')



l1 <- getLineages(reducedDim(sce_hvg.slingshot_DMSO), sce_hvg.slingshot_DMSO$GMM)
#plot_tree(pcaX, clus, l1, threeD = TRUE)
#plot_tree(reducedDim(sce_hvg.slingshot), sce_hvg.slingshot$seurat_clusters, l1, dim = 3)
plot(reducedDim(sce_hvg.slingshot_DMSO), col = brewer.pal(9,"Set1")[sce_hvg.slingshot_DMSO$seurat_clusters], asp = 1, pch = 16)
lines(l1, lwd = 3, col = 'black')

crv1 <- getCurves(l1)
crv1



## Gaussian mixture model for clustering  


cl1 <- Mclust(reducedDim(sce_hvg.slingshot_DMSO)[,1:2])$classification
cl1 <- Idents(sce_hvg.seurat_DMSO)
colData(sce_hvg.slingshot_DMSO)$GMM <- cl1

plot(reducedDim(sce_hvg.slingshot_DMSO)[,1:2], col = brewer.pal(6,"Dark2")[cl1], pch=16, asp = 1)



sce_hvg.slingshot_DMSO <- slingshot::slingshot(sce_hvg.slingshot_DMSO, clusterLabels = 'GMM', reducedDim = 'UMAP')

summary(sce_hvg.slingshot_DMSO$slingPseudotime_1)

####### Final plot to show


colors <- colorRampPalette(brewer.pal(12,'Spectral')[-6])(100)
plotcol <- colors[cut(sce_hvg.slingshot_DMSO$slingPseudotime_1, breaks=100)]

plot(reducedDims(sce_hvg.slingshot_DMSO)$UMAP, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce_hvg.slingshot_DMSO), lwd=2, col='black')


plot(reducedDims(sce_hvg.slingshot_DMSO)$UMAP, col = brewer.pal(9,'Dark2')[sce_hvg.slingshot_DMSO$seurat_clusters], pch=16, asp = 1)
lines(SlingshotDataSet(sce_hvg.slingshot_DMSO), lwd=2, type = 'lineages', col = 'black')



plot(reducedDims(sce_hvg.slingshot_DMSO)$UMAP, col = brewer.pal(9,'Dark2')[sce_hvg.slingshot_DMSO$seurat_clusters], pch=16,asp = 1)
lines(SlingshotDataSet(sce_hvg.slingshot_DMSO), lwd=2, col = 'black')
####### Final plot to show




################EPZ pseudotime



sce_hvg.seurat_EPZ_2 <- subset(sce_hvg.seurat_EPZ, idents = c("0", "1", "2", "3", "4"))
sce_hvg.slingshot_EPZ <- as.SingleCellExperiment(sce_hvg.seurat_EPZ_2)

#cl1 <- Mclust(reducedDim(sce_hvg.slingshot)[,1:2])$classification
cl1_epz <- Idents(sce_hvg.seurat_EPZ_2)
colData(sce_hvg.slingshot_EPZ)$GMM <- cl1_epz

plot(reducedDim(sce_hvg.slingshot_EPZ)[,1:2], col = brewer.pal(6,"Dark2")[cl1_epz], pch=16, asp = 1)



sce_hvg.slingshot_EPZ <- slingshot::slingshot(sce_hvg.slingshot_EPZ, clusterLabels = 'GMM', reducedDim = 'UMAP')


#saveRDS(file="sce_hvg.slingshot", sce_hvg.slingshot)
#sce_hvg.slingshot <- readRDS(file="sce_hvg.slingshot")

#summary(sce_hvg.slingshot_DMSO$slingPseudotime_1)


colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce_hvg.slingshot_EPZ$slingPseudotime_1, breaks=100)]

plot(reducedDims(sce_hvg.slingshot_EPZ)$UMAP, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce_hvg.slingshot_EPZ), lwd=2, col='black')




crv1 <- getCurves(l1)
crv1




plot(reducedDims(sce_hvg.slingshot_EPZ)$UMAP, col = brewer.pal(9,'Set1')[sce_hvg.slingshot_EPZ$seurat_clusters], pch=16, asp = 1)
lines(SlingshotDataSet(sce_hvg.slingshot_EPZ), lwd=2, type = 'lineages', col = 'black')





####### Final plot to show


colors <- colorRampPalette(brewer.pal(12,'Spectral')[-6])(100)
plotcol <- colors[cut(sce_hvg.slingshot_EPZ$slingPseudotime_1, breaks=100)]

plot(reducedDims(sce_hvg.slingshot_EPZ)$UMAP, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce_hvg.slingshot_EPZ), lwd=2, col='black')


plot(reducedDims(sce_hvg.slingshot_EPZ)$UMAP, col = brewer.pal(9,'Dark2')[sce_hvg.slingshot_EPZ$seurat_clusters], pch=16, asp = 1)
lines(SlingshotDataSet(sce_hvg.slingshot_EPZ), lwd=2, type = 'lineages', col = 'black')



plot(reducedDims(sce_hvg.slingshot_EPZ)$UMAP, col = brewer.pal(9,'Dark2')[sce_hvg.slingshot_EPZ$seurat_clusters], pch=16, asp = 1)
lines(SlingshotDataSet(sce_hvg.slingshot_EPZ), lwd=2, col = 'black')
####### Final plot to show





















