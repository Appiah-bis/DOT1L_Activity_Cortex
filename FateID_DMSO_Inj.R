
#StemID analysis
ltr <- Ltree(sc_aMI_2)
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
g <- names(ps$nodes)[ps$nodes == 19]
plotexpression(fs,y,"Nkain3",n$f,col=fcol,name="Nkain3",cluster=T,alpha=.5,types= NULL)
#types=c("DMSO24ALL","EPZ24ALL")
g

plotFateMap(y,dr,k=2,m="tsne")