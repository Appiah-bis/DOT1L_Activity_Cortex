
# Analysis of Microinjected dataset with RaceID

# For detailed explanation of each step, 
#please refer to the RaceID Manual: https://cran.r-project.org/web/packages/RaceID/vignettes/RaceID.html

# load packages.
library(RaceID)
library(FateID)
require(tsne)
require(pheatmap)
require(MASS)
require(cluster)
require(mclust)
require(flexmix)
require(lattice)
require(fpc)
require(amap)
require(RColorBrewer)
require(locfit)
require(vegan)
require(Rtsne)
require(scran)
require(DESeq2)
require(randomForest)
library(RColorBrewer)
library(stringr)
library(tidyverse)
library(reshape2)
library(ggplot2)

# Load data, initialize single cell object, perform RaceID analysis
AMI_ALL_NEW <- readRDS("AMI_ALL_NEW.RData")
sc_aMI <- SCseq(AMI_ALL_NEW)

## Gene names with correlated expressiion to any of these will be filtered out for cell type inference. 
#See: https://rdrr.io/cran/RaceID/man/filterdata.html
gm <- rownames(sc_aMI@ndata)[grep("Gm\\d+","Kcnq1","Hba-a1|Hba","Gapdh","RP23-341H2.4","A430073D23Rik", rownames(sc_aMI@ndata))]

sc_aMI <- filterdata(sc_aMI, mintotal=2500, minexpr = 10, minnumber = 5, CGenes = gm, ccor = 0.4)
fdata <- getfdata(sc_aMI)
sc_aMI <- compdist(sc_aMI,metric="spearman")
sc_aMI <- clustexp(sc_aMI,rseed = 14500)
plotsaturation(sc_aMI,disp=FALSE)
plotsaturation(sc_aMI,disp=TRUE)

plotjaccard(sc_aMI)
sc_aMI <- clustexp(sc_aMI,cln=5,sat=FALSE)
sc_aMI <- findoutliers(sc_aMI,outlg = 10)
plotbackground(sc_aMI)
plotsensitivity(sc_aMI)
plotoutlierprobs(sc_aMI)
set.seed(847)
sc_aMI_2 <- comptsne(sc_aMI)

##########
##########
# Visualize clusters | Figure 4C in manuscript, color of clusters changed to match rest of the data in Fig. 4
plotmap(sc_aMI_2,cex = 2.5)



#Extracting matrix with cellnames for Ctrl and EPZ conditions
DMSO24ALL_1<- sc_aMI@ndata[, grep('^aMIDMSO', colnames(sc_aMI@ndata))]
DMSO24ALL <- colnames(DMSO24ALL_1)
EPZ24ALL_1<- sc_aMI@ndata[, grep('^aMIEPZ', colnames(sc_aMI@ndata))]
EPZ24ALL <- colnames(EPZ24ALL_1)



# Visualize marker expression
#g <- c("Sox2", "Fabp7", "Neurog2", "Eomes", "RP23-271M19.2", "Mme", "Dpysl3", "Tubb3"), Markers to plot
#plotexpmap(sc_aMI_2, "Sox2", logsc = T, cex = 2.5) # See marker expression in all cells
plotexpmap(sc_aMI_2, "Sox2", logsc = T, cex = 2.5, cells = DMSO24ALL) # See marker expression ONLY in DMSO condition
plotexpmap(sc_aMI_2, "Sox2", logsc = T, cex = 2.5, cells = EPZ24ALL) # See marker expression ONLY in DMSO condition


# Visualze gene expression in Dotplot 

genes2plot <- c("Hes1","Hes5","Sox9","Ptprz1","Ttyh1","Fabp7","Sox2","Nusap1","Fgfr3","Hk2","Zbtb20","Nr2f1","Insm1",
                "Sox5","Eomes","Neurog2","Syt11","RP23-271M19.2","Zic1","Mme","Nkain3","A430073D23Rik","Fezf2","Pou3f2","Map2","Gap43","Dcx","Neurod6","Dpysl3","Elavl3","Tubb3","Tbr1","Bcl11b")

##########
##########
# Supplementary Figure 5C, Appiah et al. Rearranged, clusters annotated
fractDotPlot(sc_aMI, genes2plot, cluster=c(2,1,5,3,4), zsc=TRUE) # "RP23-271M19.2" is shown as Ofd1 in the manuscript


##########
##########
# Supplementary Figure 4D
plotexpmap(sc_aMI_2, "Sox2", logsc = T, cex = 2.5, cells = DMSO24ALL, n = "Sox2 | Con")
plotexpmap(sc_aMI_2, "Fabp7", logsc = T, cex = 2.5, cells = DMSO24ALL, n = "Fabp7 | Con")
plotexpmap(sc_aMI_2, "Neurog2", logsc = T, cex = 2.5, cells = DMSO24ALL, n = "Neurog2 | Con")
plotexpmap(sc_aMI_2, "Eomes", logsc = T, cex = 2.5, cells = DMSO24ALL, n = "Eomes | Con")
plotexpmap(sc_aMI_2, "RP23-271M19.2", logsc = T, cex = 2.5, cells = DMSO24ALL, n = "Ofd1 | Con")
plotexpmap(sc_aMI_2, "Mme", logsc = T, cex = 2.5, cells = DMSO24ALL, n = "Mme | Con")
plotexpmap(sc_aMI_2, "Dpysl3", logsc = T, cex = 2.5, cells = DMSO24ALL, n = "Dpysl3 | Con")
plotexpmap(sc_aMI_2, "Tubb3", logsc = T, cex = 2.5, cells = DMSO24ALL, n = "Tubb3 | Con")

plotexpmap(sc_aMI_2, "Sox2", logsc = T, cex = 2.5, cells = EPZ24ALL, n = "Sox2 | EPZ")
plotexpmap(sc_aMI_2, "Fabp7", logsc = T, cex = 2.5, cells = EPZ24ALL, n = "Fabp7 | EPZ")
plotexpmap(sc_aMI_2, "Neurog2", logsc = T, cex = 2.5, cells = EPZ24ALL, n = "Neurog2 | EPZ")
plotexpmap(sc_aMI_2, "Eomes", logsc = T, cex = 2.5, cells = EPZ24ALL, n = "Eomes | EPZ")
plotexpmap(sc_aMI_2, "RP23-271M19.2", logsc = T, cex = 2.5, cells = EPZ24ALL, n = "Ofd1 | EPZ")
plotexpmap(sc_aMI_2, "Mme", logsc = T, cex = 2.5, cells = EPZ24ALL, n = "Mme | EPZ")
plotexpmap(sc_aMI_2, "Dpysl3", logsc = T, cex = 2.5, cells = EPZ24ALL, n = "Dpysl3 | EPZ")
plotexpmap(sc_aMI_2, "Tubb3", logsc = T, cex = 2.5, cells = EPZ24ALL, n = "Tubb3 | EPZ")



