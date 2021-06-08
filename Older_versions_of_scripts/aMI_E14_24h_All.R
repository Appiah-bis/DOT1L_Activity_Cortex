sc_aMIall <- SCseq(AMI_ALL2)

gm <- rownames(sc_aMIall@ndata)[grep("Gm\\d+","Kcnq1","Hba-a1|Hba","Gapdh","RP23-341H2.4","A430073D23Rik", rownames(sc_aMIall@ndata))]



sc_aMIall <- filterdata(sc_aMIall,mintotal=2000,minexpr = 10,minnumber = 5,CGenes = gm, ccor = 0.4)
fdata <- getfdata(sc_aMIall)
sc_aMIall <- compdist(sc_aMIall,metric="logpearson")
sc_aMIall <- clustexp(sc_aMIall,rseed = 17000,clustnr = 15)
plotsaturation(sc_aMIall,disp=FALSE)
plotsaturation(sc_aMIall,disp=TRUE)

plotjaccard(sc_aMIall)
sc_aMIall <- clustexp(sc_aMIall,cln=4,sat=FALSE)
sc_aMIall <- findoutliers(sc_aMIall,outlg = 10)
plotbackground(sc_aMIall)
plotsensitivity(sc_aMIall)
plotoutlierprobs(sc_aMIall)
set.seed(245)
sc_aMIall <- comptsne(sc_aMIall)
plotmap(sc_aMIall,cex = 2.5)

##Extracting marker genes for each cluster
#marker genes for cluster 1
c1_marker <- clustdiffgenes(sc_aMIall,1,pvalue=.01)
write.csv2(c1_marker,file = "ami_ALL_c1marker.csv")

#marker genes for cluster 2
c2_marker <- clustdiffgenes(sc_aMIall,2,pvalue=.01)
write.csv2(c2_marker,file = "ami_ALL_c2marker.csv")

#marker genes for cluster 3
c3_marker <- clustdiffgenes(sc_aMIall,3,pvalue=.01)
write.csv2(c3_marker,file = "ami_ALL_c3marker.csv")

#marker genes for cluster 4
c4_marker <- clustdiffgenes(sc_aMIall,4,pvalue=.01)
write.csv2(c4_marker,file = "ami_ALL_c4marker.csv")

#Extracting matrix with cellnames for Ctrl and EPZ conditions
DMSO24ALL_1<- sc_aMIall@ndata[, grep('^aMIDMSO', colnames(sc_aMIall@ndata))]
DMSO24ALL <- colnames(DMSO24ALL_1)
EPZ24ALL_1<- sc_aMIall@ndata[, grep('^aMIEPZ', colnames(sc_aMIall@ndata))]
EPZ24ALL <- colnames(EPZ24ALL_1)
#DMSO24ALL <- c(aMIDMSO,dmso_24h_names)
#EPZ24ALL <- c(aMIEPZ,EPZ_24h_names)

#Computing contribution of each condition to proportion of cells in each cluster
clust <- sc_aMIall@cpart
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
capture.output(fisher_results,file = "fisher_results_AMIALL")

#DEG analysis - Ctrl vs EPZ
all_AMI <- diffexpnb(getfdata(sc_aMIall,n=c(DMSO24ALL,EPZ24ALL)), A=DMSO24ALL, B=EPZ24ALL)
DEGs_aMI_24h_all <- all_AMI$res[order(all_AMI$res$log2FoldChange),]











#Final_DEGS <- z$res[order(z$res$log2FoldChange),]
write.csv2(DEGs_aMI_24h_all,file = "DEGs_aMI_24h_all.csv")

write.csv2(AMI_ALL2,file = "ami_seu.csv")

s
####Differentially expressed genes in each cluster compared to all other clusters - upregulated

clustdiffgenes(EUE14ALL_2,2,pvalue=.01)

rkt <- clustdiffgenes(sc_aMIall,4,pvalue=.01)



plotexpmap2(sc_aMIall,c("Tbr1","Fezf2","Bcl11b"),n="Deep layer neurons - Fezf2/Tbr1/Ctip2",coexp = T,logsc = T,cex = 2.5,gex = "Satb2",cells = EPZ24ALL)

plotexpmap2(sc_aMIall,c("Tbr1","Fezf2","Bcl11b"),n="Deep layer neurons - Fezf2/Tbr1/Ctip2",coexp = T,logsc = T,cex = 2.5,gex = "Satb2",cells = DMSO24ALL)


plotexpmap2(sc_aMIall,c("Pax6","Sox2","Vim","Fabp7"),n="Apical Progenitors - Pax6/Sox2/Vim/Fabp7",coexp = T,logsc = T,cex = 2.5,gex = "Eomes")

plotexpmap2(sc_aMIall,c("Eomes","Insm1"),n="basal Intermediate Progenitors - Eomes/Insm1",coexp = T,logsc = T,cex = 2.5,gex = "Pax6")


library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)


