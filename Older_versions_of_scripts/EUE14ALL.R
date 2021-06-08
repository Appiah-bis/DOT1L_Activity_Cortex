#Extracting EUE DMSO and EPZ at E14
EUE14DMS <- DOT1L[which(colnames(DOT1L) %in% EUE14DMSO)]
EUE14EP <- DOT1L[which(colnames(DOT1L) %in% EUE14EPZ)]
EUE14ALL <- merge(EUE14DMS, EUE14EP, by=0)
EUE14ALL[is.na(EUE14ALL)] <- 0
rownames(EUE14ALL)=EUE14ALL$Row.names
EUE14ALL$Row.names <- NULL


EUE14ALL_2 <- SCseq(EUE14ALL)

gm <- rownames(EUE14ALL_2@ndata)[grep("Gm\\d+","Kcnq1","Hba-a1|Hba","Gapdh","RP23-341H2.4","A430073D23Rik", rownames(EUE14ALL_2@ndata))]



EUE14ALL_2 <- filterdata(EUE14ALL_2,mintotal=2500,minexpr = 10,minnumber = 5,CGenes = gm, ccor = 0.4)
fdata <- getfdata(EUE14ALL_2)
EUE14ALL_2 <- compdist(EUE14ALL_2,metric="logpearson")
EUE14ALL_2 <- clustexp(EUE14ALL_2,rseed = 17000)
plotsaturation(EUE14ALL_2,disp=FALSE)
plotsaturation(EUE14ALL_2,disp=TRUE)

plotjaccard(EUE14ALL_2)
EUE14ALL_2 <- clustexp(EUE14ALL_2,cln=4,sat=FALSE)
EUE14ALL_2 <- findoutliers(EUE14ALL_2,outlg = 10)
plotbackground(EUE14ALL_2)
plotsensitivity(EUE14ALL_2)
plotoutlierprobs(EUE14ALL_2)
set.seed(235)
EUE14ALL_2 <- comptsne(EUE14ALL_2)
plotmap(EUE14ALL_2,cex = 2.5)


##Extracting marker genes for each cluster
#marker genes for cluster 1
c1_marker <- clustdiffgenes(EUE14ALL_2,1,pvalue=.01)
write.csv2(c1_marker,file = "EUE14ALL_2_c1marker.csv")

#marker genes for cluster 2
c2_marker <- clustdiffgenes(EUE14ALL_2,2,pvalue=.01)
write.csv2(c2_marker,file = "EUE14ALL_2_c2marker.csv")

#marker genes for cluster 3
c3_marker <- clustdiffgenes(EUE14ALL_2,3,pvalue=.01)
write.csv2(c3_marker,file = "EUE14ALL_2_c3marker.csv")

#marker genes for cluster 4
c4_marker <- clustdiffgenes(EUE14ALL_2,4,pvalue=.01)
write.csv2(c4_marker,file = "EUE14ALL_2_c4marker.csv")
#marker genes for cluster 5
c5_marker <- clustdiffgenes(EUE14ALL_2,5,pvalue=.01)
write.csv2(c5_marker,file = "EUE14ALL_2_c5marker.csv")




#Computing contribution of each condition to proportion of cells in each cluster
clust <- EUE14ALL_2@cpart
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
EUE14DM <- EUE14ALL_2@ndata[, grep('^EUEE14DMSO', colnames(EUE14ALL_2@ndata))]
EUE14DM <- colnames(EUE14DM)
EUE14EP <- EUE14ALL_2@ndata[, grep('^EUEE14EPZ', colnames(EUE14ALL_2@ndata))]
EUE14EP <- colnames(EUE14EP)
##########




#Perform DEG analysis - Ctrl vs EPZ and export output as CSV
E14_EUE <- diffexpnb(getfdata(EUE14ALL_2,n=c(EUE14DM,EUE14EP)), A=EUE14DM, B=EUE14EP)
DEGs_E14_EUE_all <- E14_EUE$res[order(E14_EUE$res$log2FoldChange),]
write.csv2(DEGs_E14_EUE_all,file = "DEGs_E14_EUE_all.csv")



#Final_DEGS <- z$res[order(z$res$log2FoldChange),]
write.csv2(DEGs_aMI_24h_all,file = "DEGs_aMI_24h_all.csv")

write.csv2(AMI_ALL2,file = "ami_seu.csv")



















#DMSO24ALL <- c(aMIDMSO,dmso_24h_names)
#EPZ24ALL <- c(aMIEPZ,EPZ_24h_names)

z <- diffexpnb(getfdata(sc_1,n=c(EUE14DMSO,EUE14EPZ)), A=EUE14DMSO, B=EUE14EPZ)
head(z$res[order(z$res$padj),],80)
