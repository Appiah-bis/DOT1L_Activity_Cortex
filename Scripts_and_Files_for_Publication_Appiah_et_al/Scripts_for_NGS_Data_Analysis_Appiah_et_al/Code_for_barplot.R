# Generating barplot using targets of regulons which are in the DEGs

# Load pacckages
library(ggplot2)
library(tidyverse)


# Load DEG table
TF_gene <- read.csv2("DEGs_targeted_by_TFs.csv", sep = ",")
head(TF_gene)
p.adjust <- TF_gene$p_val_adj

##########
##########
# Plot Figure 5C, Appiah et al. 
ggplot(data = hjh, aes(x=ID, y=avg_logFC))+
  geom_bar(stat = "Identity", aes(fill = p.adjust) )+
  coord_flip()
