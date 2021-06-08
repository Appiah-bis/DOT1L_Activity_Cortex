library(tidyverse)
library(tidyr)
library(dplyr)
library(SeuratWrappers)
library(Seurat)
library(RColorBrewer)

# Heatmap of all positive markers from transit cells extracted from Electroporated cells
# Heatmap of top positive markers of AP, TC, BP, N 
#TC_ID.csv
newgene <- as.character(read.csv2("TC_ID.csv")[,1])
#tc_new <- as.character(newgene[,1])

# Subset seurat object of electroporated cells, extract expression values for only transit cells
# Select clusters to include for heatmap based on idents
# idents = 3,subset = condition == "DMSO"  idents = c("Neurons II", "Transit cells", "Intermediate Progenitors")
kty <- subset(sce_hvg.seurat) 
z <- DoHeatmap(kty,features = newgene,label = F,
               group.colors = brewer.pal(8,'Dark2')) + NoLegend()
y <- z$data %>% drop_na()
x <- y %>% group_by(Identity) %>% dplyr::select(Feature, Cell,Identity , Expression) %>%
  tidyr::spread(key = Feature, value = Expression)
#w <- y %>% dplyr::select(Feature, Cell, Expression) %>%
#  tidyr::spread(key = Cell, value = Expression) %>% column_to_rownames("Feature") %>% as.matrix()


library(ComplexHeatmap)
library(circlize)
res_list = x
type = res_list$Identity
type <- sort.int(type)
res_list <- t(as.matrix(res_list[, -c(1,2)]))
#mat_scaled = t(as.matrix(res_list[, -c(1,2)]))
#mat_scaled <- t(mat_scaled)
ha <- HeatmapAnnotation(type = type, annotation_name_side = "left",
                      show_annotation_name = TRUE,show_legend = T)

ht_list2 <- Heatmap(res_list, name = "TC Markers",
                  top_annotation = ha, 
                  show_column_names = T, row_title = NULL,show_row_dend = FALSE,
                  cluster_columns = F,show_column_dend = FALSE,row_names_side = "left"
                  ,border = T,row_names_gp = gpar(fontsize = 10),cluster_rows = T,column_gap = 5) #cluster_rows = FALSE,

draw(ht_list2, ht_gap = unit(10, "mm")) 


