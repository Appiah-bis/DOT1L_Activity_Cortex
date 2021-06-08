## install bioconductor packages
#source("https://bioconductor.org/biocLite.R")
#biocLite("scran")
#biocLite("DESeq2")
#biocLite("biomaRt")

## load required packages.
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

AMI_ALL_NEW <- readRDS("AMI_ALL_NEW.RData")
#rm(AMI_ALL_NEW)
sc_aMI <- SCseq(AMI_ALL_NEW)

gm <- rownames(sc_aMI@ndata)[grep("Gm\\d+","Kcnq1","Hba-a1|Hba","Gapdh","RP23-341H2.4","A430073D23Rik", rownames(sc_aMI@ndata))]



sc_aMI <- filterdata(sc_aMI,mintotal=2500,minexpr = 10,minnumber = 5,CGenes = gm, ccor = 0.4)
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
sc_aMI_3 <- compfr(sc_aMI)
sc_aMI <- compumap(sc_aMI)
plotmap(sc_aMI,cex = 2.5,um=T)

#Extracting matrix with cellnames for Ctrl and EPZ conditions
DMSO24ALL_1<- sc_aMI@ndata[, grep('^aMIDMSO', colnames(sc_aMI@ndata))]
DMSO24ALL <- colnames(DMSO24ALL_1)
EPZ24ALL_1<- sc_aMI@ndata[, grep('^aMIEPZ', colnames(sc_aMI@ndata))]
EPZ24ALL <- colnames(EPZ24ALL_1)

#Check coexpression of markers
coexpmk <- c("RP23-271M19.2","Syt11","A230046K03Rik","Ccdc36","Mme")


genes2 <- c("Hes1","Hes5","Sox9","Ptprz1","Ttyh1","Fabp7","Sox2","Nusap1","Fgfr3","Hk2","Zbtb20","Nr2f1","Insm1",
            "Sox5","Eomes","Neurog2","Syt11","RP23-271M19.2","Zic1","Mme","Nkain3","A430073D23Rik","Fezf2","Pou3f2","Map2","Gap43","Dcx","Neurod6","Dpysl3","Elavl3","Tubb3","Tbr1","Bcl11b")
genes3 <- c("Hes1","Hes5","Sox9","Ptprz1","Ttyh1","Fabp7","Sox2","Add1","Add2","Add3","Zbtb20","Nr2f1","Insm1",
                      "Sox5","Eomes","Neurog2","Syt11","RP23-271M19.2","Zic1","Dcx","Neurod6","Dpysl3","Elavl3","Tubb3","Tbr1","Bcl11b")
# Visualze gene expression with Dotplots

fractDotPlot(sc_aMI, genes3, cluster=c(2,1,5,3,4), zsc=TRUE)




 #StemID analysis
 ltr <- Ltree(sc_aMI_2)
 ltr <- compentropy(ltr)
 ltr <- projcells(ltr,cthr=5,nmode=F,fr=F)
 ltr <- projback(ltr,pdishuf=100)
 ltr <- lineagegraph(ltr)
 ltr <- comppvalue(ltr,pthr=0.1)
 plotgraph(ltr,scthr=0.2,showCells=FALSE)
 x <- compscore(ltr,scthr=0.2)
 plotdistanceratio(ltr)
 plotspantree(ltr,cex = 2)
 plotprojections(ltr)

 #Inspecting pseudotemporal gene expression changes
 n <- cellsfromtree(ltr,c(2,1,5,3,4))

 
 x <- getfdata(ltr@sc)
 library(FateID)
 dr  <- compdr(x, z=NULL, m=c("tsne","cmd","dm","lle","umap"), k=c(2,3), lle.n=30, dm.sigma="local", dm.distance="euclidean", tsne.perplexity=30, seed=12345)
 tar <- c(2,1,5,3,4)

 
 fs  <- filterset(x,n=n$f)
 s1d <- getsom(fs,nb=1000,alpha=.5)
 ps  <- procsom(s1d,corthr=.85,minsom=3)
 y    <- ltr@sc@cpart[n$f]
 fcol <- ltr@sc@fcol
 plotheatmap(ps$nodes.z,xpart=y,xcol=fcol,ypart=unique(ps$nodes),xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
 plotheatmap(ps$all.z,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
 plotheatmap(ps$all.e,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
 plotheatmap(ps$all.b,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
 g <- names(ps$nodes)[ps$nodes == 1]
 plotexpression(fs,y,"Sox9",n$f,col=fcol,name="Sox9",cluster=T,alpha=.5,types=c("DMSO24ALL","EPZ24ALL"))
g
plotFateMap(y,fs,k=2,m="tsne")

tyz <- c(DMSO24ALL,EPZ24ALL)


clust <- sc_aMI_2@cpart
clust_pts <- data.frame(names=names(clust),cluster=clust)
clust_pts$condition <- substr(clust_pts$names, start = 0, stop = 6)
table_clust_pts <- table(clust_pts$condition, clust_pts$cluster)
table_clust_pts_1 <- table_clust_pts[1:2,]

#Fishers Exact test for contribution of condition (Ctrl/EPZ) to proportions of cells in each cluster

datalist <- list()
data24h <- data.frame("CellsIn"=c(1,1), "CellsOut"=c(1,1))

for (j in 1:dim(table_clust_pts_1)[2]) {
  for (i in 1:dim(table_clust_pts_1)[1]) {
    data24h$CellsIn[i] <- table_clust_pts_1[i,j]
    data24h$CellsOut[i] <- (sum(table_clust_pts_1[i,]) - table_clust_pts_1[i,j])
    datalist[[j]] <- data24h
    rownames(datalist[[j]]) <- rownames(table_clust_pts_1)
  }
}

fisher_results <- list()

for (i in 1:length(datalist)) {
  fisher_results[[i]] <- fisher.test(datalist[[i]])
  
}
#Export results from Fishers Exact test
#capture.output(fisher_results,file = "fisher_results_AMI_ALL")







