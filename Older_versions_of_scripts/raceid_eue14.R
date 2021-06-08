#Extracting EUE DMSO and EPZ at E14
EUE14DMS <- DOT1L[which(colnames(DOT1L) %in% EUE14DMSO)]
EUE14EP <- DOT1L[which(colnames(DOT1L) %in% EUE14EPZ)]
EUE14ALL <- merge(EUE14DMS, EUE14EP, by=0)
EUE14ALL[is.na(EUE14ALL)] <- 0
rownames(EUE14ALL)=EUE14ALL$Row.names
EUE14ALL$Row.names <- NULL


#data3 <- readRDS("EUE14ALL.RData")
sc_EUE14 <- SCseq(data3)

gm <- rownames(sc_EUE14@ndata)[grep("Gm\\d+","Kcnq1","Hba-a1|Hba","Gapdh","RP23-341H2.4","A430073D23Rik", rownames(sc_EUE14@ndata))]



sc_EUE14 <- filterdata(sc_EUE14,mintotal=2500,CGenes = gm, ccor = 0.4)
fdata <- getfdata(sc_EUE14)
sc_EUE14 <- compdist(sc_EUE14,metric="spearman")
sc_EUE14 <- clustexp(sc_EUE14,rseed = 30000)
plotsaturation(sc_EUE14,disp=FALSE)
plotsaturation(sc_EUE14,disp=TRUE)

plotjaccard(sc_EUE14)
sc_EUE14 <- clustexp(sc_EUE14,cln=4,sat=FALSE)
sc_EUE14 <- findoutliers(sc_EUE14,outlg = 20)
plotbackground(sc_EUE14)
plotsensitivity(sc_EUE14)
plotoutlierprobs(sc_EUE14)
set.seed(5000)
sc_EUE14 <- comptsne(sc_EUE14)
plotmap(sc_EUE14,cex = 2.5)


##Extracting marker genes for each cluster
#marker genes for cluster 1
c1_marker <- clustdiffgenes(sc_EUE14,1,pvalue=.01)
write.csv2(c1_marker,file = "EUE14ALL_2_c1marker.csv")

#marker genes for cluster 2
c2_marker <- clustdiffgenes(sc_EUE14,2,pvalue=.01)
write.csv2(c2_marker,file = "EUE14ALL_2_c2marker.csv")

#marker genes for cluster 3
c3_marker <- clustdiffgenes(sc_EUE14,3,pvalue=.01)
write.csv2(c3_marker,file = "EUE14ALL_2_c3marker.csv")

#marker genes for cluster 4
c4_marker <- clustdiffgenes(sc_EUE14,4,pvalue=.01)
write.csv2(c4_marker,file = "EUE14ALL_2_c4marker.csv")
#marker genes for cluster 5
c5_marker <- clustdiffgenes(sc_EUE14,5,pvalue=.01)
write.csv2(c5_marker,file = "EUE14ALL_2_c5marker.csv")




#Computing contribution of each condition to proportion of cells in each cluster
clust <- sc_EUE14@cpart
clust_pts <- data.frame(names=names(clust),cluster=clust)
clust_pts$condition <- substr(clust_pts$names, start = 0, stop = 9)
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
capture.output(fisher_results,file = "fisher_results_EUE14ALL_2")

##########
EUE14DM <- sc_EUE14@ndata[, grep('^EUEE14DMSO', colnames(sc_EUE14@ndata))]
EUE14DM <- colnames(EUE14DM)
EUE14EP <- sc_EUE14@ndata[, grep('^EUEE14EPZ', colnames(sc_EUE14@ndata))]
EUE14EP <- colnames(EUE14EP)
##########




#Perform DEG analysis - Ctrl vs EPZ and export output as CSV
E14_EUE <- diffexpnb(getfdata(sc_EUE14,n=c(EUE14DM,EUE14EP)), A=EUE14DM, B=EUE14EP)
DEGs_E14_EUE_all <- E14_EUE$res[order(E14_EUE$res$log2FoldChange),]
write.csv2(DEGs_E14_EUE_all,file = "DEGs_E14_EUE_all.csv")



#Final_DEGS <- z$res[order(z$res$log2FoldChange),]
write.csv2(DEGs_aMI_24h_all,file = "DEGs_aMI_24h_all.csv")

write.csv2(AMI_ALL2,file = "ami_seu.csv")



















#DMSO24ALL <- c(aMIDMSO,dmso_24h_names)
#EPZ24ALL <- c(aMIEPZ,EPZ_24h_names)

z <- diffexpnb(getfdata(sc_1,n=c(EUE14DMSO,EUE14EPZ)), A=EUE14DMSO, B=EUE14EPZ)
head(z$res[order(z$res$padj),],80)