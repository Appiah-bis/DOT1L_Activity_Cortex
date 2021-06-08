#Extract cells from Ctrl and treatment condition for StemID
DMSO24ALL_1<- sc_aMI_2@expdata[, grep('^aMIDMSO', colnames(sc_aMI_2@expdata))]

DMSO24ALL_1 <- SCseq(DMSO24ALL_1)

gm <- rownames(DMSO24ALL_1@ndata)[grep("Gm\\d+","Kcnq1","Hba-a1|Hba","Gapdh","RP23-341H2.4","A430073D23Rik", rownames(DMSO24ALL_1@ndata))]



DMSO24ALL_1 <- filterdata(DMSO24ALL_1,mintotal=2500,minexpr = 10,minnumber = 5,CGenes = gm, ccor = 0.4)
fdata <- getfdata(DMSO24ALL_1)
DMSO24ALL_1 <- compdist(DMSO24ALL_1,metric="spearman")
DMSO24ALL_1 <- clustexp(DMSO24ALL_1,rseed = 14500)
plotsaturation(DMSO24ALL_1,disp=FALSE)
plotsaturation(DMSO24ALL_1,disp=TRUE)

plotjaccard(DMSO24ALL_1)
DMSO24ALL_1 <- clustexp(DMSO24ALL_1,cln=5,sat=FALSE)
DMSO24ALL_1 <- findoutliers(DMSO24ALL_1,outlg = 10)
plotbackground(DMSO24ALL_1)
plotsensitivity(DMSO24ALL_1)
plotoutlierprobs(DMSO24ALL_1)
set.seed(847)
DMSO24ALL_2 <- comptsne(DMSO24ALL_1)
DMSO24ALL_1 <- compumap(DMSO24ALL_1)
plotmap(DMSO24ALL_1,cex = 2.5,um=T)



#StemID analysis
ltr <- Ltree(DMSO24ALL_2)
ltr <- compentropy(ltr)
ltr <- projcells(ltr,cthr=5,nmode=F,fr=F)
ltr <- projback(ltr,pdishuf=100)
ltr <- lineagegraph(ltr)
ltr <- comppvalue(ltr,pthr=0.1)
plotgraph(ltr,scthr=0.2,showCells=FALSE,showMap = T,cex = 2.5)
x <- compscore(ltr,scthr=0.2)
plotdistanceratio(ltr)
plotspantree(ltr,cex = 2.5)#,projections = T
plotprojections(ltr)

#Inspecting pseudotemporal gene expression changes
n <- cellsfromtree(ltr,c(2,1,5,3,4))


x <- getfdata(ltr@sc)
library(FateID)




#tar <- c(2,5,4)


fs  <- filterset(x,n=n$f)
s1d <- getsom(fs,nb=1000,alpha=.5)
ps  <- procsom(s1d,corthr=.85,minsom=3)
y    <- ltr@sc@cpart[n$f]
fcol <- ltr@sc@fcol
plotheatmap(ps$nodes.z,xpart=y,xcol=fcol,ypart=unique(ps$nodes),xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
plotheatmap(ps$all.z,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
plotheatmap(ps$all.e,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
plotheatmap(ps$all.b,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
g <- names(ps$nodes)[ps$nodes == 18]
plotexpression(fs,y,"Sox2",n$f,col=fcol,name="Sox2",cluster=T,alpha=.5,types= NULL)
#types=c("DMSO24ALL","EPZ24ALL")
g

plotFateMap(y,dr,k=2,m="tsne")
