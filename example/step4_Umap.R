#----
#library packages
library(pacman)
pacman::p_load(ggplot2,dplyr,tidyr,Matrix, Seurat, tidyverse)

setwd("/your_payh/re") #change

# The following naming corresponds to:
# All - all cell
# Center - core cell
# Normal - intermediate cell
# Border - marginal cell

scRNAlist2 <- list()
scRNA <- readRDS("scRNAanno.rds") ## Annotated RDS object list of four populations: all cell, core cell, intermadiate cell, marginal cell
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
DefaultAssay(scRNA_CNB) <- "integrated"
scRNA_CNB <- ScaleData(scRNA_CNB,
                       #vars.to.regress = c("S.Score", "G2M.Score","percent_mito"),
                       features =rownames(scRNA_CNB))
scRNA_CNB <- RunPCA(scRNA_CNB)
scRNA_CNB <- FindNeighbors(scRNA_CNB, reduction = "pca")
scRNA_CNB <- RunUMAP(scRNA_CNB,dims=1:30)
scRNA_CNB <- RunTSNE(scRNA_CNB,dims=1:30)

#(2)integration
#scRNAlist_reintegrate
scRNAlist <- SplitObject(scRNA, split.by = "CellClass")
for(x in 1:length(scRNAlist)){
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


#UMAP
umap_theme <- theme( 
  axis.line=element_blank(), 
  axis.text.x=element_blank(), 
  axis.text.y=element_blank(), 
  axis.ticks=element_blank(), 
  axis.title.x=element_blank(), 
  axis.title.y=element_blank(), 
  panel.background=element_blank(), 
  panel.border=element_blank(), 
  panel.grid.major=element_blank(), 
  panel.grid.minor=element_blank()
) 

my_cols <- c("activated Hepatic stellate cell" = "#DB7093", 
             "resting Hepatic stellate cell" = "#8B864E",
             "Central vein endothelial cell" = "#C1FFC1",
             "Hepatic artery endothelial cell" = "#C0937E",
             "Liver sinusoidal endothelial cell" = "#FFD700",
             "Portal vein endothelial cell" = "#A4D38E",
             "Lymphatic endothelial cell" = "#1C1C1C",
             "Hepatocyte" = "#FF7F24",
             "Cholangiocyte" = "#CD0000",
             "cNK" = "#3477A9",
             "lrNK" = "#BDA7CB",#liver resident NK
             "Erythrocyte" = "#684797",
             "B cell" =  "#F19294",
             #"Plasma cell" = "#C0FF3E",
             "Plasma cell" = "#8B008B",
             "cDC1"= "#4A9D47",
             "cDC2" = "#8B5A2B",
             "pDC" = "#006400",
             "Neutrophil progenitor" = "#FFA500",
             "Mature neutrophil" = "#FFDAB9",
             "Monocyte" = "#696969",
             "Kupffer cell" = "#BBFFFF",
             "MAIT" = "#B9D3EE",
             "Mast cell" = "#00CED1",
             "Naive T cell" = "#123456",
             "CD8+ T" = "#6495ED",
             "NKT" = "#EE82EE",
             "Proliferative cell" = "#8B0000",
             "Unknown" = "#ece9e9")

my_cols2 <- my_cols[order(as.integer(names(my_cols)))]
scales::show_col(my_cols2)

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


