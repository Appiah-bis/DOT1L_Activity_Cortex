#Extract cells from Ctrl and treatment condition for StemID
EPZ24ALL_1<- sc_aMI_2@expdata[, grep('^aMIEPZ', colnames(sc_aMI_2@expdata))]

EPZ24ALL_1 <- SCseq(EPZ24ALL_1)

gm <- rownames(EPZ24ALL_1@ndata)[grep("Gm\\d+","Kcnq1","Hba-a1|Hba","Gapdh","RP23-341H2.4","A430073D23Rik", rownames(EPZ24ALL_1@ndata))]



EPZ24ALL_1 <- filterdata(EPZ24ALL_1,mintotal=2500,minexpr = 10,minnumber = 5,CGenes = gm, ccor = 0.4)
fdata <- getfdata(EPZ24ALL_1)
EPZ24ALL_1 <- compdist(EPZ24ALL_1,metric="spearman")
EPZ24ALL_1 <- clustexp(EPZ24ALL_1,rseed = 14500)
plotsaturation(EPZ24ALL_1,disp=FALSE)
plotsaturation(EPZ24ALL_1,disp=TRUE)

plotjaccard(EPZ24ALL_1)
EPZ24ALL_1 <- clustexp(EPZ24ALL_1,cln=5,sat=FALSE)
EPZ24ALL_1 <- findoutliers(EPZ24ALL_1,outlg = 10)
plotbackground(EPZ24ALL_1)
plotsensitivity(EPZ24ALL_1)
plotoutlierprobs(EPZ24ALL_1)
set.seed(847)
EPZ24ALL_2 <- comptsne(EPZ24ALL_1)
EPZ24ALL_1 <- compumap(EPZ24ALL_1)
plotmap(EPZ24ALL_1,cex = 2.5,um=T)



#StemID analysis
ltr <- Ltree(EPZ24ALL_2)
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
n <- cellsfromtree(ltr,c(4,5,2,3,1))


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
g <- names(ps$nodes)[ps$nodes == 4]
plotexpression(fs,y,"Sox2",n$f,col=fcol,name="Sox2",cluster=T,alpha=.5,types= NULL)
#types=c("DMSO24ALL","EPZ24ALL")
g

plotFateMap(y,dr,k=2,m="tsne")