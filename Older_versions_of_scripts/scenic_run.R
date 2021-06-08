

#Load required packages
#suppressPackageStartupMessages({
 # library(SCENIC)
  #library(AUCell)
  #library(RcisTarget)
  #library(SCopeLoomR)
#})

#options(width=200)

#Load pre-processed expression matrix -- containing 800 highly variable genes from Seurat analysis
#sce_hvg_1 <- readRDS("sce_hvg.RDS")
#exprMat <- counts(sce_hvg_1)
#cellInfo <- colData(sce_hvg_1)

### Initialize settings
#library(SCENIC)


#org <- "mgi" # or hgnc, or dmel
#dbDir <- "~/" # RcisTarget databases location
#myDatasetTitle <- "MyData" # choose a name for your analysis
#data(defaultDbNames)
#dbs <- defaultDbNames[[org]]
#scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=5) 
#scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.RDS"
#saveRDS(scenicOptions, file="int/scenicOptions.RDS") 

### Co-expression network
#genesKept <- geneFiltering(exprMat, scenicOptions)


#genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
 #                          minCountsPerGene=3*.01*ncol(exprMat),
  #                         minSamples=ncol(exprMat)*.01)
#exprMat_filtered <- exprMat[genesKept, ]
#runCorrelation(exprMat_filtered, scenicOptions)
#exprMat_filtered_log <- log2(exprMat_filtered+1) 
#runGenie3(exprMat_filtered_log, scenicOptions)

### Build and score the GRN
#exprMat_log <- log2(exprMat+1)
#scenicOptions <- readRDS("int/scenicOptions.Rds")
#scenicOptions@settings$verbose <- TRUE
#scenicOptions@settings$nCores <- 6
#scenicOptions@settings$seed <- 123
#scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] 

#runSCENIC_1_coexNetwork2modules(scenicOptions)
#runSCENIC_2_createRegulons(scenicOptions,coexMethod = "top50perTarget") # Toy run settings
#runSCENIC_3_scoreCells(scenicOptions, exprMat_log)


library(KernSmooth)
library(RColorBrewer)
dens2d <- bkde2D(tSNE_scenic$Y, 1, gridsize = c(dim(tSNE_scenic$Y)[1],dim(tSNE_scenic$Y)[1]))$fhat
image(dens2d, col=brewer.pal(9, "YlOrBr"), axes=FALSE)
contour(dens2d, add=TRUE, nlevels=5, drawlabels=FALSE)

dens2d_b_new1 <- dens2d_b[-indice_cond1,-indice_cond1]
dens2d_b_new2 <- dens2d_b[-indice_cond2,-indice_cond2]
dens2d_b_new3 <- dens2d_b[-indice_cond3,-indice_cond3]


tSNE_scenic_new <- tSNE_scenic$Y[-indice_cond1,]
indice_cond1 <- which(!rownames(dens2d_b) %in% dmso_24h_names)
indice_cond2 <- which(!rownames(dens2d_b) %in% dmso_48h_names)
indice_cond3 <- which(!rownames(dens2d_b) %in% EPZ_24h_names)


# append cell names from SCENIC tSNE as rownames and colnames of dens2d and rename
dens2d_b <- dens2d
cell_names <- rownames(tSNE_scenic$Y)
rownames(dens2d_b) <- cell_names
colnames(dens2d_b) <- cell_names
dms24_GRN <- rownames(dens2d_b)%in%dmso_2









