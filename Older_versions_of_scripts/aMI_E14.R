sc_ami <- SCseq(AMI_ALL)

gm <- rownames(sc_ami@ndata)[grep("Gm\\d+","Kcnq1","Hba-a1|Hba","Gapdh","RP23-341H2.4","A430073D23Rik", rownames(sc_ami@ndata))]



sc_ami <- filterdata(sc_ami,mintotal=2000,minexpr = 10,minnumber = 5,CGenes = gm, ccor = 0.4)
fdata <- getfdata(sc_ami)
sc_ami <- compdist(sc_ami,metric="logpearson")
sc_ami <- clustexp(sc_ami,rseed = 17000,clustnr = 15)
plotsaturation(sc_ami,disp=FALSE)
plotsaturation(sc_ami,disp=TRUE)

plotjaccard(sc_ami)
sc_ami <- clustexp(sc_ami,cln=3,sat=FALSE)
sc_ami <- findoutliers(sc_ami,outlg = 10)
#sc_ami <- findoutliers(sc_ami)
plotbackground(sc_ami)
plotsensitivity(sc_ami)
plotoutlierprobs(sc_ami)
set.seed(272)
sc_ami <- comptsne(sc_ami)
plotmap(sc_ami,cex = 2.5)

#Extracting cell names from aMI cells
aMIDMSOE14<- sc_ami@ndata[, grep('^aMIDMSOE14', colnames(sc_ami@ndata))]
aMIDMSO <- colnames(aMIDMSOE14)
aMIEPZE14<- sc_ami@ndata[, grep('^aMIEPZE14', colnames(sc_ami@ndata))]
aMIEPZ <- colnames(aMIEPZE14)

z <- diffexpnb(getfdata(sc_ami,n=c(aMIDMSO,aMIEPZ)), A=aMIDMSO, B=aMIEPZ)
ami_DEGnew_padj <- z$res[order(z$res$padj),]
DEG <- head(z$res[order(z$res$padj),],100)

##Extracting marker genes for each cluster
#marker genes for cluster 1
c1_marker <- clustdiffgenes(sc_ami,1,pvalue=.01)
write.csv2(c1_marker,file = "ami_exp2_c1marker.csv")

#marker genes for cluster 2
c2_marker <- clustdiffgenes(sc_ami,2,pvalue=.01)
write.csv2(c2_marker,file = "ami_exp2_c2marker.csv")

#marker genes for cluster 3
c3_marker <- clustdiffgenes(sc_ami,3,pvalue=.01)
write.csv2(c3_marker,file = "ami_exp2_c3marker.csv")

#marker genes for cluster 4
c4_marker <- clustdiffgenes(sc_ami,4,pvalue=.01)
write.csv2(c4_marker,file = "ami_exp2_c4marker.csv")








