##  Slingshot for pseudotime analysis 
# Load packages


library(BiocStyle)
library(readr)
library(data.table)
library(dplyr)
library(tools)
library(SingleCellExperiment)
library(HDF5Array)
library(mvoutlier)
library(scran)
library(scater)
library(SC3)
library(slingshot)
library(stringr)
library(org.Mm.eg.db)
library(mclust)
library(RColorBrewer)
library(Seurat)
library(umap)
library(latexpdf)
library(slingshot)


###################################################################
# Slingshot pseudotime analysis for Control (DMSO) cells only #####
###################################################################

###### load RData files for analysis: sce_hvg.seurat_DMSO_2.RData
sce_hvg.seurat <- readRDS("sce_hvg.seurat_DMSO_2.RData")
sce_hvg.slingshot_DMSO <- as.SingleCellExperiment(sce_hvg.seurat)
cl1 <- Mclust(reducedDim(sce_hvg.slingshot_DMSO)[,1:2])$classification
cl1 <- Idents(sce_hvg.seurat_DMSO)
colData(sce_hvg.slingshot_DMSO)$GMM <- cl1

plot(reducedDim(sce_hvg.slingshot_DMSO)[,1:2], col = brewer.pal(6,"Dark2")[cl1], pch=16, asp = 1)


sce_hvg.slingshot_DMSO <- slingshot::slingshot(sce_hvg.slingshot_DMSO, clusterLabels = 'GMM', reducedDim = 'UMAP')

summary(sce_hvg.slingshot_DMSO$slingPseudotime_1)

####### Plot result


# plotcol <- colors[cut(sce_hvg.slingshot_DMSO$slingPseudotime_1, breaks=100)]
# 
# plot(reducedDims(sce_hvg.slingshot_DMSO)$UMAP, col = plotcol, pch=16, asp = 1)
# lines(SlingshotDataSet(sce_hvg.slingshot_DMSO), lwd=2, col='black')


plot(reducedDims(sce_hvg.slingshot_DMSO)$UMAP, col = brewer.pal(9,'Dark2')[sce_hvg.slingshot_DMSO$seurat_clusters], pch=16, asp = 1)
lines(SlingshotDataSet(sce_hvg.slingshot_DMSO), lwd=2, type = 'lineages', col = 'black')



plot(reducedDims(sce_hvg.slingshot_DMSO)$UMAP, col = brewer.pal(9,'Dark2')[sce_hvg.slingshot_DMSO$seurat_clusters], pch=16,asp = 1)
lines(SlingshotDataSet(sce_hvg.slingshot_DMSO), lwd=2, col = 'black')
####### 




###################################################################
# Slingshot pseudotime analysis for EPZ cells only                #
###################################################################

###### load RData files for analysis: sce_hvg.seurat_DMSO_2.RData

sce_hvg.seurat_EPZ <- readRDS("sce_hvg.seurat_EPZ_2.RData")
sce_hvg.slingshot_EPZ <- as.SingleCellExperiment(sce_hvg.seurat_EPZ)
cl1 <- Mclust(reducedDim(sce_hvg.slingshot_EPZ)[,1:2])$classification
cl1 <- Idents(sce_hvg.seurat_EPZ)
colData(sce_hvg.slingshot_EPZ)$GMM <- cl1

plot(reducedDim(sce_hvg.slingshot_EPZ)[,1:2], col = brewer.pal(6,"Dark2")[cl1], pch=16, asp = 1)


sce_hvg.slingshot_EPZ <- slingshot::slingshot(sce_hvg.slingshot_EPZ, clusterLabels = 'GMM', reducedDim = 'UMAP')

summary(sce_hvg.slingshot_EPZ$slingPseudotime_1)

# plotcol <- colors[cut(sce_hvg.slingshot_EPZ$slingPseudotime_1, breaks=100)]
# plot(reducedDims(sce_hvg.slingshot_EPZO)$UMAP, col = plotcol, pch=16, asp = 1)
# lines(SlingshotDataSet(sce_hvg.slingshot_EPZ), lwd=2, col='black')

plot(reducedDims(sce_hvg.slingshot_EPZ)$UMAP, col = brewer.pal(6,'Dark2')[sce_hvg.slingshot_EPZ$seurat_clusters], pch=16, asp = 1)
lines(SlingshotDataSet(sce_hvg.slingshot_EPZ), lwd=2, type = 'lineages', col = 'black')



plot(reducedDims(sce_hvg.slingshot_EPZ)$UMAP, col = brewer.pal(6,'Dark2')[sce_hvg.slingshot_EPZ$seurat_clusters], pch=16,asp = 1)
lines(SlingshotDataSet(sce_hvg.slingshot_EPZ), lwd=2, col = 'black')
