# Load packages
library(DOSE)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)


# GO analysis using genes downregulated in TTS compared to BP only in DMSO condition

# Load data file named Decreased_In_C3_comp_IP_DMSO_electroporated.csv for this analysis
symbols_1 <- as.character(read.csv2("Decreased_In_C3_comp_IP_DMSO_electroporated.csv")[,1]) #TTS_vs_AP_DEG_UP.csv
gene <- mapIds(org.Mm.eg.db, symbols_1, 'ENTREZID', 'SYMBOL')

ego_1 <- enrichGO(gene          = gene,
                  OrgDb         = org.Mm.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  readable      = TRUE)

dotplot(ego_1, showCategory=20, orderBy="GeneRatio",title = "TTS vs BP (DMSO) - Decreased Genes")

############################################################


# GO analysis using genes upregulated in TTS compared to BP only in DMSO condition
# Load data file named Increased_In_C3_comp_IP_DMSO_electroporated.csv for this analysis
symbols_1 <- as.character(read.csv2("Increased_In_C3_comp_IP_DMSO_electroporated.csv")[,1]) #TTS_vs_AP_DEG_UP.csv
gene <- mapIds(org.Mm.eg.db, symbols_1, 'ENTREZID', 'SYMBOL')

ego_1 <- enrichGO(gene          = gene,
                  OrgDb         = org.Mm.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  readable      = TRUE)

dotplot(ego_1, showCategory=20, orderBy="GeneRatio",title = "TTS vs BP (DMSO) - Increased Genes")


###########################################################


# GO analysis using genes downregulated in TTS compared to Neurons only in DMSO condition

# Load data file named Decreased_In_C3_comp_Neu_DMSO_electroporated.csv for this analysis
symbols_1 <- as.character(read.csv2("Decreased_In_C3_comp_Neu_DMSO_electroporated.csv")[,1]) #TTS_vs_AP_DEG_UP.csv
gene <- mapIds(org.Mm.eg.db, symbols_1, 'ENTREZID', 'SYMBOL')

ego_1 <- enrichGO(gene          = gene,
                  OrgDb         = org.Mm.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  readable      = TRUE)


dotplot(ego_1, showCategory=20,orderBy="GeneRatio",title = "TTS vs Neurons (DMSO) - Decreased Genes")
##########################

# GO analysis using genes upregulated in TTS compared to Neurons only in DMSO condition
# Load data file named Increased_In_C3_comp_Neu_DMSO_electroporated.csv for this analysis
symbols_1 <- as.character(read.csv2("Increased_In_C3_comp_Neu_DMSO_electroporated.csv")[,1]) #TTS_vs_AP_DEG_UP.csv
gene <- mapIds(org.Mm.eg.db, symbols_1, 'ENTREZID', 'SYMBOL')

ego_1 <- enrichGO(gene          = gene,
                  OrgDb         = org.Mm.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  readable      = TRUE)


dotplot(ego_1, showCategory=20,orderBy="GeneRatio",title = "TTS vs Neurons (DMSO) - Increased Genes")
##########################

