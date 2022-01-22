# Analysis of Electroporation datset using single cell gene regulatory network inference and clustering (SCENIC)
# Perform seurat analysis and then export the seurat objects for control and treated cells for  this analysis

# Initial prepaprations
# All required objects and data files are provided in this repository. Use  the file names provided here to copy them to your working directory
######
# You may load previous analysis where I used 10 PCs for clustering electroporated data. The neurons cluster is split into two (NI and NII).
# Only TF activites for NII is shown as it covers both neuronal populations here
#Load the objects: "sce_hvg.seurat_DMSO_new"   "sce_hvg.seurat_EPZ_new" 
# Load the environment: "SCENIC_RUN" 
#Install all packages and their dependencies
# IMPORTANT: All results will go to your active working directory, see folder "output"
#########

########

# For details about the method, see -  Aibar et al. (2017) SCENIC: single-cell regulatory network inference and clustering. Nature Methods. doi:
# DOI: https://doi.org/10.1038/nmeth.4463
# Or visit: https://github.com/aertslab/SCENIC
##########

# Load packages
library(GEOquery)
library(data.table)
library(Biobase)
library(SingleCellExperiment)
library(SCENIC)
library(GENIE3)
library(AUCell)
library(RcisTarget)
library(SCopeLoomR)
library(Seurat)
library(pheatmap)
library(DT)
library(pheatmap)
library(R2HTML)
library(doMC)
library(doRNG)




# perform seurat analysis of electroporation data, extract counts for cells in each condition
#see markdown file: Code_for_Electroporation_Dataset.Rmd
# Alternatively, load seurat objects provided
#sce_hvg_1 <- load(sce_hvg.seurat_DMSO_new)
#sce_hvg_2 <- load(sce_hvg.seurat_EPZ_new)
sce_hvg_1 <- subset(sce_hvg.seurat_DMSO_new, subset= condition=="DMSO")
sce_hvg_2 <- subset(sce_hvg.seurat_EPZ_new, subset= condition=="EPZ")
exprMat <-as.matrix(GetAssayData(sce_hvg_1, slot = "counts"))
cellInfo <- WhichCells(sce_hvg_1)


### Initialize settings
library(SCENIC)


org <- "mgi" # or hgnc, or dmel
dbDir <- "~/" # RcisTarget databases location
myDatasetTitle <- "MyData" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=5) 



scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.RDS"
saveRDS(scenicOptions, file="int/scenicOptions.RDS") 

### Co-expression network
#genesKept <- geneFiltering(exprMat, scenicOptions)

genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered_log, scenicOptions)

### Build and score the GRN
exprMat_log <- log2(exprMat+1)
#exprMat_filtered <- exprMat[genesKept, ]
#rm(exprMat)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- readRDS("int/scenicOptions.Rds")
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 6
scenicOptions@settings$seed <- 123
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] 

runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions,coexMethod = "top50perTarget") # Toy run settings
runSCENIC_3_scoreCells(scenicOptions, exprMat_log)

# Export:
scenicOptions@fileNames$output["loomFile",] <- "myAnalysis.loom"
export2scope(scenicOptions, myAnalysis)


#Binarize activity?
aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_log)
savedSelections <- shiny::runApp(aucellApp)
newThresholds <- savedSelections$thresholds
scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
saveRDS(scenicOptions, file="int/scenicOptions.Rds")
runSCENIC_4_aucell_binarize(scenicOptions)

### Exploring output 
# Check files in folder 'output'
# .loom file @ http://scope.aertslab.org

# output/Step2_MotifEnrichment_preview.html in detail/subset:
motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")
tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="Rad21"]
viewMotifs(tableSubset) 

# output/Step2_regulonTargetsInfo.tsv in detail: 
regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
tableSubset <- regulonTargetsInfo[TF=="Stat6" & highConfAnnot==TRUE]
viewMotifs(tableSubset) 


library(KernSmooth)
library(RColorBrewer)
dens2d <- bkde2D(tSNE_scenic$Y, 1, gridsize = c(dim(tSNE_scenic$Y)[1],dim(tSNE_scenic$Y)[1]))$fhat
image(dens2d_b, col=brewer.pal(9, "YlOrBr"), axes=FALSE)
contour(dens2d_b, add=TRUE, nlevels=5, drawlabels=FALSE)


cellInfo2 <- data.frame(seuratCluster=Idents(sce_hvg_1))


regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCluster <- sapply(split(rownames(cellInfo), cellInfo$seuratCluster),
                                    function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCluster_Scaled <- t(scale(t(regulonActivity_byCluster), center = T, scale=T))

pheatmap::pheatmap(regulonActivity_byCluster_Scaled, #fontsize_row=3, 
                   color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-3, 3, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color=NA)
# filename="regulonActivity_byCellType.pdf", width=10, height=20)

topRegulators <- reshape2::melt(regulonActivity_byCluster_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
viewTable(topRegulators)











