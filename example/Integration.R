library(pacman)
pacman::p_load(ggplot2,dplyr,tidyr,Matrix, gtools, ggsignif,ggsci,
               stringr,Seurat,cowplot, clustree, ggpubr, scales,
               tibble,patchwork,SCpubr,ROGUE, Hmisc, rstatix,
               Scillus,magrittr,tidyverse,reticulate, plyr,
               ggsci,reshape2,scDblFinder,SingleR,celldex,cluster, kBET)

setwd("/public/home/wutong/Project/PJ3_HumanLiver_scRNAseq_Databse/Result/Figure2/F2_V3/") #change
use_condaenv("/public/home/wutong/Software/miniconda3/envs/R4.1.1/bin/python")
py_module_available("sklearn")

#----
#== Datasets Integrate ==#

#== 
scRNA <- readRDS("/public/home/wutong/Project/PJ3_HumanLiver_scRNAseq_Databse/Result/step4_annotation/scRNAanno.rds")
Cellclass <- scRNA@meta.data$CellClass
names(Cellclass) <- rownames(scRNA@meta.data)
Idents(scRNA) <- Cellclass
scRNAlist <- list()
scRNAlist[[1]] <- subset(x = scRNA,subset = (CellClass == "Center"))
scRNAlist[[2]] <- subset(x = scRNA,subset = (CellClass == "Normal"))
scRNAlist[[3]] <- subset(x = scRNA,subset = (CellClass == "Border"))



#----
#== uamp ==#
