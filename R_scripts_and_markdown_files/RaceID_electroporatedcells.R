library(RaceID)
require(RColorBrewer)



DOT1L <- readRDS("Bismark2.RData")
sc_1 <- SCseq(DOT1L)
rm(DOT1L)
gm <- rownames(sc_1@ndata)[grep("Gm\\d+", rownames(sc_1@ndata))]

sc_1 <- filterdata(sc_1,mintotal=2500,CGenes = gm, ccor = 0.4)
fdata <- getfdata(sc_1)
sc_1 <- compdist(sc_1,metric="logpearson")
sc_1 <- clustexp(sc_1,rseed = 17000)
plotsaturation(sc_1,disp=FALSE)
plotsaturation(sc_1,disp=TRUE)

plotjaccard(sc_1)
sc_1 <- clustexp(sc_1,cln=4,sat=FALSE)
sc_1 <- findoutliers(sc_1)
plotbackground(sc_1)
plotsensitivity(sc_1)
plotoutlierprobs(sc_1)
set.seed(347)
sc_1 <- comptsne(sc_1)
plotmap(sc_1,cex = 1.5)



#Extracting cell names from EUE E14 cells

EUE_14_DMSO<- sc_1@ndata[, grep('^EUEE14DMSO', colnames(sc_1@ndata))]
EUE14DMSO <- colnames(EUE_14_DSMSO)
EUE_14_EPZ<- sc_1@ndata[, grep('^EUEE14EPZ', colnames(sc_1@ndata))]
EUE14EPZ <- colnames(EUE_14_EPZ)

##########
EUE_12_ALL <- DOT1L[, grep('^EUEE12', colnames(DOT1L))]
EUE12_ALL <- colnames(EUE_12_ALL)


##########
EUE_14_DSMSO_2 <- DOT1L[, grep('^EUEE14DMSO', colnames(DOT1L))]
EUE14DMSO_2 <- colnames(EUE_14_DSMSO_2)
EUE_14_EPZ_2 <- DOT1L[, grep('^EUEE14EPZ', colnames(DOT1L))]
EUE14EPZ_2 <- colnames(EUE_14_EPZ_2)
##########


###### E16
EUE_16_DSMSO<- sc_1@ndata[, grep('^EUEE16DMSO', colnames(sc_1@ndata))]
EUE16DMSO <- colnames(EUE_16_DSMSO)
EUE_16_EPZ<- sc_1@ndata[, grep('^E16EPZ', colnames(sc_1@ndata))]
EUE16EPZ <- colnames(EUE_16_EPZ)

#######E12
EUE_12_DSMSO<- sc_1@ndata[, grep('^EUEE12DMSO', colnames(sc_1@ndata))]
EUE12DMSO <- colnames(EUE_12_DSMSO)



#####Creating expression matrix for only micrinjected cells 
#AMI_DMSO_E14<- sc_1@expdata[, grep('^aMIDMSOE14', colnames(sc_1@expdata))]
#AMI_EPZ_E14 <- sc_1@expdata[, grep('^aMIEPZE14', colnames(sc_1@expdata))]
#rm(AMI_DMSO_E14)
#rm(AMI_EPZ_E14)

AMI_DMSO14 <- DOT1L[, grep('^aMIDMSOE14', colnames(DOT1L))]
AMI_EPZ14 <- DOT1L[, grep('^aMIEPZE14', colnames(DOT1L))]

AMI_ALL <- merge(AMI_DMSO14, AMI_EPZ14, by=0)
AMI_ALL[is.na(AMI_ALL)] <- 0
rownames(AMI_ALL)=AMI_ALL$Row.names
AMI_ALL$Row.names <- NULL

#Extracting DMSO and EPZ conditions from first injection experiment
inj_1_dmso24 <- dot1l[which(colnames(dot1l) %in% dmso_24h_names)]
colnames(inj_1_dmso24) <- paste("aMIDMSO", colnames(inj_1_dmso24), sep = "_")

inj_1_epz24 <- dot1l[which(colnames(dot1l) %in% EPZ_24h_names)]
colnames(inj_1_epz24) <- paste("aMIEPZ", colnames(inj_1_epz24), sep = "_")

inj_1_all <- merge(inj_1_dmso24, inj_1_epz24, by=0)
inj_1_all[is.na(inj_1_all)] <- 0
rownames(inj_1_all)=inj_1_all$Row.names
inj_1_all$Row.names <- NULL

#colnames(inj_1_epz24) <- paste("Sub", colnames(m2), sep = "_")


#Merging DMSO and EPZ 24h from all microinjection experiments

AMI_ALL2 <- merge(inj_1_all, AMI_ALL, by=0)
AMI_ALL2[is.na(AMI_ALL2)] <- 0
rownames(AMI_ALL2)=AMI_ALL2$Row.names
AMI_ALL2$Row.names <- NULL

#Extracting EUE DMSO and EPZ at E14
EUE14DMS <- DOT1L[which(colnames(DOT1L) %in% EUE14DMSO)]
EUE14EP <- DOT1L[which(colnames(DOT1L) %in% EUE14EPZ)]
EUE14ALL <- merge(EUE14DMS, EUE14EP, by=0)
EUE14ALL[is.na(EUE14ALL)] <- 0
rownames(EUE14ALL)=EUE14ALL$Row.names
EUE14ALL$Row.names <- NULL









#AMI_ALL





#k <- rownames(AMI_EPZ14)

z <- diffexpnb(getfdata(sc_1,n=c(EUE14DMSO,EUE14EPZ)), A=EUE14DMSO, B=EUE14EPZ)
head(z$res[order(z$res$log2FoldChange),],50)
real.degs <- head(z$res[order(z$res$padj),],100)

