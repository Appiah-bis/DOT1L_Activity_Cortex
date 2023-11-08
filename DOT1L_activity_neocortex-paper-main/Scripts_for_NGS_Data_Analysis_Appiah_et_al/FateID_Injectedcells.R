#Extract cells from Ctrl and treatment condition for StemID
library(RaceID)
library(RColorBrewer)

# load preprocessed RData object: sc_aMI_2
sc_aMI_2 = readRDS("sc_AMI.RData")
# subset DMSO treated cells
DMSO_aMI<- sc_aMI_2@expdata[, grep('^aMIDMSO', colnames(sc_aMI_2@expdata))]
DMSO_aMI <- SCseq(DMSO_aMI)

gm <- rownames(DMSO_aMI@ndata)[grep("Gm\\d+","Kcnq1","Hba-a1|Hba","Gapdh","RP23-341H2.4","A430073D23Rik", rownames(DMSO_aMI@ndata))]



DMSO_aMI <- filterdata(DMSO_aMI,mintotal=2500,minexpr = 10,minnumber = 5,CGenes = gm, ccor = 0.4)
fdata <- getfdata(DMSO_aMI)
DMSO_aMI <- compdist(DMSO_aMI,metric="spearman")
DMSO_aMI <- clustexp(DMSO_aMI,rseed = 14500)
plotsaturation(DMSO_aMI,disp=FALSE)
plotsaturation(DMSO_aMI,disp=TRUE)

plotjaccard(DMSO_aMI)
DMSO_aMI <- clustexp(DMSO_aMI,cln=5,sat=FALSE)
DMSO_aMI <- findoutliers(DMSO_aMI,outlg = 10)
plotbackground(DMSO_aMI)
plotsensitivity(DMSO_aMI)
plotoutlierprobs(DMSO_aMI)
set.seed(847)
DMSO_aMI_1<- comptsne(DMSO_aMI)
DMSO_aMI_2 <- compumap(DMSO_aMI)
plotmap(DMSO_aMI_1,cex = 2.5)
plotmap(DMSO_aMI_2,cex = 2.5,um=T)


# genes2 <- c("Hes1","Hes5","Sox9","Ptprz1","Ttyh1","Fabp7","Sox2","Nusap1","Fgfr3","Hk2","Zbtb20","Nr2f1","Insm1",
#             "Sox5","Eomes","Neurog2","Fezf2","Pou3f2","Map2","Gap43","Dcx","Neurod6","Dpysl3","Elavl3","Tubb3","Tbr1","Bcl11b")
# 
# # Visualze gene expression with Dotplots
# 
# fractDotPlot(DMSO_aMI_1, genes2, cluster=c(2,1,5,3,4), zsc=TRUE)

#StemID analysis
ltr <- Ltree(DMSO_aMI_1)
ltr <- compentropy(ltr)
ltr <- projcells(ltr,cthr = 5,nmode=F,fr=F)
ltr <- projback(ltr,pdishuf=100)
ltr <- lineagegraph(ltr)
ltr <- comppvalue(ltr, pthr = 0.1)
plotgraph(ltr,scthr=0.2,showCells=FALSE,showMap = T,cex = 2.5)
x <- compscore(ltr,scthr=0.2)
plotdistanceratio(ltr)
plotspantree(ltr,cex = 2.5)#,projections = T


#Inspecting pseudotemporal gene expression changes
n <- cellsfromtree(ltr,c(2,1,5,3,4))

x <- getfdata(ltr@sc)

library(FateID)

fs  <- filterset(x,n=n$f)
s1d <- getsom(fs,nb=1000,alpha=.5)
ps  <- procsom(s1d,corthr=.85,minsom=3)
y    <- ltr@sc@cpart[n$f]
fcol <- ltr@sc@fcol
plotheatmap(ps$nodes.z,xpart=y,xcol=fcol,ypart=unique(ps$nodes),xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
plotheatmap(ps$all.z,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
plotheatmap(ps$all.e,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
plotheatmap(ps$all.b,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
g <- names(ps$nodes)[ps$nodes == 4]

# Visualize pseudotemporal expression of selected markers
plotexpression(fs,y,"Sox2",n$f,col=fcol,name="Sox2",cluster=T,alpha=.5,types= NULL)
plotexpression(fs,y,"RP23-341H2.4",n$f,col=fcol,name="Ofd1",cluster=T,alpha=.5,types= NULL)
plotexpression(fs,y,"Eomes",n$f,col=fcol,name="Eomes",cluster=T,alpha=.5,types= NULL)
plotexpression(fs,y,"Tubb3",n$f,col=fcol,name="Tubb3",cluster=T,alpha=.5,types= NULL)

