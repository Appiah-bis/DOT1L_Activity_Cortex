ndata <- as.matrix(sc_1@ndata * min(sc_1@counts)) + 0.1
ndata <- data.frame(ndata)
cpart <- sc_1@cpart
names(cpart) <- colnames(ndata)
uniq <- sort(unique(cpart))
head(uniq)
my_list <- list()
for (i in 1:length(uniq)){
  x <- names(cpart[cpart == uniq[i]])
  my_list[[i]] <- x #assign(paste("clust", i, sep = "_"),x)
  names(my_list)[[i]] <- paste("clust", uniq[i], sep="_")
}
my_list2 <- list()
for (i in 1:length(my_list)){
  #cat("running",names(my_list)[i],"\n")
  xd <- diffexpnb(ndata, my_list[[i]], method="per-condition", colnames(ndata[,!(colnames(ndata) %in% my_list[[i]])]), norm=F,vfit=sc@background$vfit)
  xd <- xd$res
  xd <- xd[order(xd$pval),]
  #assign(paste("diff_clust", i, sep="_"), xd)
  my_list2[[i]] <- xd
  names(my_list2)[[i]] <- paste("clust", uniq[i], "vs_all", sep="_")
}
my_list3 <- list()
for (i in 1:length(my_list2)){
  my_list3[[i]] <- my_list2[[i]][my_list2[[i]]$log2FoldChange < 0,]
  names(my_list3)[[i]] <- names(my_list2)[[i]]
}
my_list7 <- list()
for (i in 1:length(my_list2)){
  my_list7[[i]] <- my_list2[[i]][my_list2[[i]]$log2FoldChange > 0,]
  names(my_list7)[[i]] <- names(my_list2)[[i]]
}
#### Julia macrophages per-condition diffexpnb
Rpl <- rownames(sc_1@ndata)[grep("Rpl", rownames(sc_1@ndata))]
Rps <- rownames(sc_1@ndata)[grep("Rps", rownames(sc_1@ndata))]
Gm <- rownames(sc_1@ndata)[grep("Gm\\d+", rownames(sc_1@ndata))]
Rp <- rownames(sc_1@ndata)[grep("RP\\d+", rownames(sc_1@ndata))]
Rik <- rownames(sc_1@ndata)[grep("Rik", rownames(sc_1@ndata))]
out <- c(Rpl, Rps, Gm, Rik, Rp)
my_list4 <- list()
for (i in 1:length(my_list3)) {
  my_list4[[i]] <- my_list3[[i]][!(rownames(my_list3[[i]]) %in% out),]
  names(my_list4)[[i]] <- names(my_list3)[[i]]
}
my_list5 <- lapply(my_list3, function(x){
  x[order(x$log2FoldChange),]
}
)
my_list6 <- lapply(my_list4, function(x){
  x[order(x$log2FoldChange),]
}
)