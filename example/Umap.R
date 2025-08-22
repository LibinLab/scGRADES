#----
#library packages
library(pacman)
pacman::p_load(ggplot2,dplyr,tidyr,Matrix, Seurat, tidyverse)

setwd("/public/home/wutong/Project/PJ3_HumanLiver_scRNAseq_Databse/Result/Figure2/F2_V3/") #change

scRNAlist2 <- list()
scRNA <- readRDS("/public/home/wutong/Project/PJ3_HumanLiver_scRNAseq_Databse/Result/step4_annotation/scRNAanno.rds")
scRNAlist_tmp <- SplitObject(scRNA, split.by = "CellClass")
scRNAtmp <- scRNAlist_tmp[[2]]
scRNAlist <- SplitObject(scRNAtmp, split.by = "samplebatch")

for(i in 1:length(scRNAlist)){
#scRNAlist[[i]] <- NormalizeData(scRNAlist[[i]])
scRNAlist[[i]] <- FindVariableFeatures(scRNAlist[[i]], selection.method = "vst",nfeatures = 5000)
}

features <- SelectIntegrationFeatures(object.list = scRNAlist, nfeatures = 2000)
scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist,dims = 1:30)     
scRNA_CNB <- IntegrateData(anchorset = scRNA.anchors, k.weight = 25)
scRNA_CNB <- ScaleData(scRNA_CNB,
                       #vars.to.regress = c("S.Score", "G2M.Score","percent_mito"),
                       features =rownames(scRNA_CNB))
scRNA_CNB <- RunPCA(scRNA_CNB)
#scRNA_CNB <- RunPCA(scRNA_CNB, features=VariableFeatures(scRNA_CNB),verbose = FALSE)
scRNA_CNB <- FindNeighbors(scRNA_CNB, reduction = "pca")
#scRNA_CNB <- FindClusters(scRNA_CNB, resolution = i)
scRNA_CNB <- RunUMAP(scRNA_CNB,dims=1:30)
scRNA_CNB <- RunTSNE(scRNA_CNB,dims=1:30)
# setwd("/public/home/wutong/Project/PJ3_HumanLiver_scRNAseq_Databse/Result/Figure2/F2_V2/")
# saveRDS(scRNA_CNB, "P1_2_scRNAlist_anno_UMAP_border.rds")

#(2)分割降维聚类
#scRNAlist_reintegrate
scRNAlist <- SplitObject(scRNA, split.by = "CellClass")
for(x in 2:length(scRNAlist)){
  scRNAlist[[x]] <- RunPCA(scRNAlist[[x]], features=VariableFeatures(scRNAlist[[x]]),verbose = FALSE)
  scRNAlist[[x]] <- FindNeighbors(scRNAlist[[x]], reduction = "pca")
  scRNAlist[[x]] <- RunUMAP(scRNAlist[[x]],dims=1:50)
  scRNAlist[[x]] <- RunTSNE(scRNAlist[[x]],dims=1:50,check_duplicates = FALSE)
}

scRNAlist[[1]]@meta.data$CellType <- factor(scRNAlist[[1]]@meta.data$CellType, 
                                  levels = c("CD4+ T","CD8+ T", "Treg","NKT","cNK","lrNK","B cell",
                                             "pDC", "cDC1","cDC2", "Monocyte", "Kupffer cell","Mast cell",
                                             "Plasma cell","Mature neutrophil","Immature neutrophil",
                                             "Fibroblast","Stellate cell","Hepatocyte","Cholangiocyte",
                                             "Erythrocyte","Liver sinusoidal endothelial cell", 
                                             "Hepatic artery endothelial cell",
                                             "Portal vein endothelial cell",
                                             "Proliferating cell","Unknown"))
scRNAlist[[2]]@meta.data$CellType <- factor(scRNAlist[[2]]@meta.data$CellType, 
                                           levels = c("CD4+ T","CD8+ T", "Treg","NKT","cNK","lrNK","B cell",
                                                      "pDC", "cDC1","cDC2", "Monocyte", "Kupffer cell","Mast cell",
                                                      "Plasma cell","Mature neutrophil","Immature neutrophil",
                                                      "Fibroblast","Stellate cell","Hepatocyte","Cholangiocyte",
                                                      "Erythrocyte","Liver sinusoidal endothelial cell", 
                                                      "Hepatic artery endothelial cell",
                                                      "Portal vein endothelial cell",
                                                      "Proliferating cell","Unknown"))
scRNAlist[[3]]@meta.data$CellType <- factor(scRNAlist[[3]]@meta.data$CellType, 
                                           levels = c("CD4+ T","CD8+ T", "Treg","NKT","cNK","lrNK","B cell",
                                                      "pDC", "cDC1","cDC2", "Monocyte", "Kupffer cell","Mast cell",
                                                      "Plasma cell","Mature neutrophil","Immature neutrophil",
                                                      "Fibroblast","Stellate cell","Hepatocyte","Cholangiocyte",
                                                      "Erythrocyte","Liver sinusoidal endothelial cell", 
                                                      "Hepatic artery endothelial cell",
                                                      "Portal vein endothelial cell",
                                                      "Proliferating cell","Unknown"))
scRNAlist[[4]]@meta.data$CellType <- factor(scRNAlist[4]]@meta.data$CellType, 
                                            levels = c("CD4+ T","CD8+ T", "Treg","NKT","cNK","lrNK","B cell",
                                                       "pDC", "cDC1","cDC2", "Monocyte", "Kupffer cell","Mast cell",
                                                       "Plasma cell","Mature neutrophil","Immature neutrophil",
                                                       "Fibroblast","Stellate cell","Hepatocyte","Cholangiocyte",
                                                       "Erythrocyte","Liver sinusoidal endothelial cell", 
                                                       "Hepatic artery endothelial cell",
                                                       "Portal vein endothelial cell",
                                                       "Proliferating cell","Unknown"))


colors <- c("#DB7093","#FFDEAD","#C1FFC1","#C0937E",
            "#FFD700","#A4D38E","#1C1C1C","#FF7F24",
            "#CD0000","#3477A9","#BDA7CB","#684797",
            "#FFFF00","#C0FF3E","#4A9D47","#8B5A2B",
            "#FFA500","#BBFFFF","#B9D3EE","#00CED1",
            "#696969","#F19294","#6495ED","#EE82EE",
            "#8B0000","#CFCFCF")

P1_C_recluster <- DimPlot(scRNAlist[[1]], 
                group.by = "CellType",
                label = F, 
                cols = colors,
                pt.size = 0.00001,
                repel = T,
                raster=FALSE) +
  labs(title= "") +
  umap_theme +
  theme(panel.border = element_rect(color="grey",fill=NA)) 

P1_N_recluster <- DimPlot(scRNAlist[[2]], 
                group.by = "CellType",
                label = F, 
                cols =  colors,
                pt.size = 0.00001,
                repel = T,
                raster=FALSE) +
  labs(title= "") +
  umap_theme +
  theme(panel.border = element_rect(color="grey",fill=NA)) 
P1_B_recluster <- DimPlot(scRNAlist[[3]], 
                group.by = "CellType",
                label = F, 
                cols = colors,
                pt.size = 0.00001,
                repel = T,
                raster=FALSE) +
  labs(title= "") +
  umap_theme +
  theme(panel.border = element_rect(color="grey",fill=NA)) 
P1_All <- DimPlot(scRNAlist[[4]], 
                          group.by = "CellType",
                          label = F, 
                          cols = colors,
                          pt.size = 0.00001,
                          repel = T,
                          raster=FALSE) +
  labs(title= "") +
  umap_theme +
  theme(panel.border = element_rect(color="grey",fill=NA)) 

P1_umap_int_recluster <- P1_All + P1_C_recluster + P1_N_recluster+ P1_B_recluster
ggsave("P1_3_umap_int_re_cluster_V2.pdf", plot = P1_umap_int_recluster , width = 20, height = 13)
ggsave("P1_3_umap_int_re_cluster_large_V2.pdf", plot = P1_umap_int_recluster , width = 33, height = 26,limitsize = FALSE,dpi=6000)
saveRDS(scRNAlist, "P1_scRNAlist_ACNB.rds")


