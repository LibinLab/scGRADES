#----
#library packages
library(pacman)
pacman::p_load(ggplot2,dplyr,tidyr,Matrix, gtools, ggsignif,ggsci,
               stringr,Seurat,cowplot, clustree, ggpubr, scales,
               tibble,patchwork,SCpubr,ROGUE, Hmisc, rstatix,
               Scillus,magrittr,tidyverse,reticulate, plyr,
               ggsci,reshape2,scDblFinder,SingleR,celldex,cluster, kBET)

setwd("/your_path/re/") #change
use_condaenv("/yourpath/miniconda3/envs/R4.1.1/bin/python")
py_module_available("sklearn")

#----
#==down sample==#
# Load the integrated scRNA object for running scGRADES, 
# with metadata containing population information for three groups
scRNA <- readRDS("scRNA.rds")
Cellclass <- scRNA@meta.data$CellClass
names(Cellclass) <- rownames(scRNA@meta.data)
Idents(scRNA) <- Cellclass
scRNAlist <- list()
scRNAlist2 <- list()

# Perform downsampling on three cell types, split the data, and conduct clustering
scRNAlist[[1]] <- subset(x = scRNA,subset = (CellPopulation == "core"))
scRNAlist[[2]] <- subset(x = scRNA,subset = (CellPopulation == "intermediate"))
scRNAlist[[3]] <- subset(x = scRNA,subset = (CellPopulation == "marginal"))

scRNAlist2[[1]] <- subset(scRNA, downsample = 20487) 
scRNAlist2[[2]] <- subset(scRNAlist[[1]], downsample = 20487)
scRNAlist2[[3]] <- subset(scRNAlist[[2]], downsample = 20487)
scRNAlist2[[4]] <- subset(scRNAlist[[3]], downsample = 20487)


scRNAlist <- list()
scRNAlist <- scRNAlist2
scRNAlist[[1]]@meta.data <- scRNAlist[[1]]@meta.data %>% 
  dplyr::rename(CellPopulation2 = CellPopulation) %>% mutate(CellPopulation = "All")

SStable <- list()
smp_meta <- list()
Smp_meta <- list()
for(x in 1:length(scRNAlist)){
  for(i in seq(0.1,1.5,0.1)){
    n = i * 10
    DefaultAssay(scRNAlist[[x]]) <- "integrated" #用重整合数据需要指定assay
    scRNAlist[[x]] <- RunPCA(scRNAlist[[x]], features=VariableFeatures(scRNAlist[[x]]),verbose = FALSE)
    scRNAlist[[x]] <- FindNeighbors(scRNAlist[[x]], reduction = "pca")
    scRNAlist[[x]] <- FindClusters(scRNAlist[[x]], resolution = i)
    smp_meta[[n]] <-  scRNAlist[[x]]@meta.data %>%
      tidygraph::select(str_c("integrated_snn_res.", i))
  }
  Smp_meta[[x]] <- smp_meta
  #scRNAlist[[x]] <- RunUMAP(scRNAlist[[x]],dims=1:50)
}
saveRDS(scRNAlist,"Res15_subset.rds")
saveRDS(Smp_meta,"Res15_CNB_meta.rds")

SStable <- list()
smp_meta <- list()
Smp_meta <- list()
for(x in 1:length(scRNAlist)) {
   for (i in seq(0.1, 1.5, 0.1)) {
     n = i * 10
     smp_meta[[n]] <-  scRNAlist[[x]]@meta.data %>%
       tidygraph::select(str_c("integrated_snn_res.", i))
     Smp_meta[[x]] <- smp_meta
   }
 }



for(y in 1:length(scRNAlist)){
  name <- head(scRNAlist[[y]]@meta.data$CellPopulation,1)
  tmpSmpmeta <- data.matrix(data.frame(Smp_meta[[y]]))
  assign(name,tmpSmpmeta)
}

#----
#== CNB NMI ==# 
##Cell class NMI
repl_python()
import sklearn
from sklearn import metrics
from sklearn.metrics.cluster import normalized_mutual_info_score

tmpNMI_all = list()
tmpNMI_core = list()
tmpNMI_intermediate = list()
tmpNMI_marginal = list()
NMIlist_all = list()
NMIlist_core = list()
NMIlist_intermediate = list()
NMIlist_marginal = list() #change
for x, y in zip(range(0,14), range(1,15)):
  tmpNMI_all = normalized_mutual_info_score(r.all[:,x], r.all[:,y])
  tmpNMI_core = normalized_mutual_info_score(r.core[:,x], r.core[:,y]) #change
  tmpNMI_intermediate = normalized_mutual_info_score(r.intermediate[:,x], r.intermediate[:,y])
  tmpNMI_marginal = normalized_mutual_info_score(r.marginal[:,x], r.marginal[:,y])
  NMIlist_all.append(tmpNMI_all)
  NMIlist_center.append(tmpNMI_core) #change
  NMIlist_normal.append(tmpNMI_intermediate)
  NMIlist_marginal.append(tmpNMI_marginal)
quit

NMI_list <- list()
NMI_list$all <- data.frame(py$NMIlist_all) 
NMI_list$core <- data.frame(py$NMIlist_core) 
NMI_list$intermediate <- data.frame(py$NMIlist_intermediate)
NMI_list$marginal <- data.frame(py$NMIlist_marginal)

##Plot
NMI_smp_list2 <- list()
WeightNMI_smp_list <- list()

for(z in 1:length(NMI_list)){ 
  name <- head(scRNAlist[[z]]@meta.data$CellPopulation,1)
  tmp_NMI_smp <- NMI_list[[z]] %>%
    set_names(str_c(seq(0.1,1.4,0.1),"_",seq(0.2,1.5,0.1))) %>%
    t() %>%
    data.frame() %>%
    mutate(Lable = rep(paste(name,"",sep = ""),ncol(NMI_list[[z]]))) 
  #colnames(tmp_NMI_smp) <- NULL
  tmp_NMI_smp[,3] <- rownames(tmp_NMI_smp)
  rownames(tmp_NMI_smp) <- NULL
  names(tmp_NMI_smp) <- c("NMI","CellPopulation","Resolution")
  NMI_smp_list2[[z]] <- tmp_NMI_smp
}
AllNMI <- do.call(rbind,NMI_smp_list2)

AllNMI$CellClass <- as.factor(AllNMI$CellClass)
AllNMI$CellClass <- factor(AllNMI$CellClass, levels = c("all","core","intermediate","marginal"))

ptmp <- ggplot(data = AllNMI,
               mapping = aes(x = Resolution, y = NMI, group= CellClass, color = CellClass)) +
  theme_bw() +
  geom_line() +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(panel.grid = element_blank()) +
  ylim(0.6, 1) +
  theme(axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10)) +
  scale_fill_manual(values = c("#F3C861","#e98307","#91BF6E","#b8c8db")) +
  scale_color_manual(values = c("#F3C861","#e98307","#91BF6E","#b8c8db")) +
  #scale_fill_manual(values = c("#3675b6","#72a2d5","#bcbdbf")) +
  #scale_color_manual(values = c("#3675b6","#72a2d5","#bcbdbf")) +
  labs(y = "NMI", x = "Resolution") 

ggsave("P2_CNB_NMI_0.6.pdf", plot = ptmp, width = 3.5, height = 2.7) #change
write.table(AllNMI,
            file = "P2_CNB_NMI.txt", #change
            quote =F,
            sep = "\t")

#----
#== CNB homogeneity score ==#
repl_python()
import sklearn
from sklearn.metrics.cluster import homogeneity_score
tmpC_all = list()
tmpC_center = list()
tmpC_normal = list()
tmpC_border = list()
Clist_all = list()
Clist_center = list() #change
Clist_normal = list()
Clist_border = list()
for x, y in zip(range(0,14), range(1,15)):
  tmpC_all = homogeneity_score(r.All[:,x], r.All[:,y])
  tmpC_center = homogeneity_score(r.Center[:,x], r.Center[:,y])#change
  tmpC_normal = homogeneity_score(r.Normal[:,x], r.Normal[:,y])
  tmpC_border = homogeneity_score(r.Border[:,x], r.Border[:,y])
  Clist_all.append(tmpC_all)
  Clist_center.append(tmpC_center)
  Clist_normal.append(tmpC_normal)
  Clist_border.append(tmpC_border)#change
quit

C_list <- list()
C_list$All <- data.frame(py$Clist_all) 
C_list$Center <- data.frame(py$Clist_center) #change
C_list$Normal <- data.frame(py$Clist_normal)
C_list$Border <- data.frame(py$Clist_border)
#...... write all sample 

C_smp_list2 <- list()
for(z in 1:length(C_list)){ #change
  name <- head(scRNAlist[[z]]@meta.data$CellClass,1)
  tmp_C_smp <- C_list[[z]] %>%
    set_names(str_c(seq(0.1,1.4,0.1),"_",seq(0.2,1.5,0.1))) %>%
    t() %>%
    data.frame() %>%
    mutate(Lable = rep(paste("",name,sep = ""),ncol(C_list[[z]]))) 
  tmp_C_smp[,3] <- rownames(tmp_C_smp)
  rownames(tmp_C_smp) <- NULL
  names(tmp_C_smp) <- c("Homogeneity_score","CellClass","Resolution")
  C_smp_list2[[z]] <- tmp_C_smp
}
AllC <- do.call(rbind, C_smp_list2)
AllC$CellClass <- as.factor(AllC$CellClass)
AllC$CellClass <- factor(AllC$CellClass, levels = c("All","Center","Normal","Border"))
##PLOT
ptmp_C <- ggplot(data = AllC,
               mapping = aes(x = Resolution, y = Homogeneity_score, group= CellClass, color = CellClass)) +
  theme_bw() +
  geom_line() +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(panel.grid = element_blank()) +
  ylim(0.8, 1.01) +
  theme(axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10)) +
  scale_fill_manual(values = c("#F3C861","#e98307","#91BF6E","#b8c8db")) +
  scale_color_manual(values = c("#F3C861","#e98307","#91BF6E","#b8c8db")) +
  labs(y = "Homogeneity score", x = "Resolution") 

ggsave("P3_CNB_HS_0.8.pdf", plot = ptmp_C, width = 3.5, height = 2.7) #change
write.table(AllC,
            file = "P3_CNB_HS.txt", #change
            quote =F,
            sep = "\t")

#----
#== Cell Type ASW ==#

# The following naming corresponds to:
# All - all cell
# Center - core cell
# Normal - intermediate cell
# Border - marginal cell

library(pacman)
pacman::p_load(ggplot2,tidyverse,Hmisc,
               gtools,ggsignif,ggsci,
               ggpubr,reshape2,rstatix,scales)

for (n in 1:length(scRNAlist)) {
  PCA <- Embeddings(object = scRNAlist[[n]], reduction = "pca")
  meta <- list()
  for(i in seq(0.1,1.5,0.1)){
    m = i * 10
    DefaultAssay(scRNAlist[[n]]) <- "integrated"
    #scRNAlist[[n]] <- FindNeighbors(scRNAlist[[n]], reduction = "pca")
    #scRNAlist[[n]] <- FindClusters(scRNAlist[[n]], resolution = i)
    meta[[m]] <-  scRNAlist[[n]]@meta.data %>% 
      select(str_c("integrated_snn_res.", i))
  }
  metaM <- data.matrix(data.frame(meta))
  name <- head(scRNAlist[[n]]@meta.data$CellClass,1)
  name1 <- paste("DataMetalist_",name,sep = "")
  name2 <- paste("DataPCAlist_",name,sep = "")
  assign(name1,metaM)
  assign(name2,PCA)
}

repl_python()
import sklearn
from sklearn import metrics
SSlist_All = list()
SSlist_Center = list()
SSlist_Normal = list()
SSlist_Border = list()
for n in range(15):
  tmpSS_All = metrics.silhouette_score(r.DataPCAlist_All, r.DataMetalist_All[:,n])
  tmpSS_Center = metrics.silhouette_score(r.DataPCAlist_Center, r.DataMetalist_Center[:,n])
  tmpSS_Normal = metrics.silhouette_score(r.DataPCAlist_Normal, r.DataMetalist_Normal[:,n])
  tmpSS_Border = metrics.silhouette_score(r.DataPCAlist_Border, r.DataMetalist_Border[:,n])
  SSlist_All.append(tmpSS_All)
  SSlist_Center.append(tmpSS_Center)
  SSlist_Normal.append(tmpSS_Normal)
  SSlist_Border.append(tmpSS_Border)
quit

ss_list <- list()
ss_list$All <- data.frame(py$SSlist_All)
ss_list$Center <- data.frame(py$SSlist_Center) #change
ss_list$Normal <- data.frame(py$SSlist_Normal)
ss_list$Border <- data.frame(py$SSlist_Border)

sslist2 <- list()
meansslist <- list()
medsslist <- list()
meanssctlist <- list()
medssctlist <- list()
for(z in 1:length(ss_list)){
  name <- head(scRNAlist[[z]]@meta.data$CellClass,1)
  tmp_ss <- ss_list[[z]] %>%
    set_names(str_c("res", seq(0.1,1.5,0.1), sep = "")) %>%
    t() %>%
    data.frame() %>%
    mutate(Dataset = rep(name)) 
  tmp_ss <- tmp_ss %>% mutate(Resolution = row.names(tmp_ss))
  rownames(tmp_ss) <- NULL
  colnames(tmp_ss) <- c("ASW","CellType","Resolution") 
  tmp_ss <- tmp_ss %>% mutate(CellTypeASW = (ASW + 1) / 2)
  mean_ss <- mean(tmp_ss$ASW) 
  med_ss <- median(tmp_ss$ASW)
  mean_ctss <- mean(tmp_ss$CellTypeASW) 
  med_ctss <- median(tmp_ss$CellTypeASW)
  sslist2[[z]] <- tmp_ss
  meansslist[[z]] <- mean_ss
  medsslist[[z]] <- med_ss
  meanssctlist[[z]] <- mean_ctss
  medssctlist[[z]] <- med_ctss
}
Allss <- do.call(rbind,sslist2)
Meanss <- do.call(rbind,meansslist)
Medss <- do.call(rbind,medsslist)
Meanssct <- do.call(rbind,meanssctlist)
Medssct <- do.call(rbind,medssctlist)

write.table(Allss,
            file = "P4_CellTypeASW.txt", #change
            quote =F,
            sep = "\t")
write.table(Meanss,
            file = "P4_CellTypeASW_mean.txt", #change
            quote =F,
            sep = "\t")
write.table(Medss,
            file = "P4_CellTypeASW_med.txt", #change
            quote =F,
            sep = "\t")
write.table(Meanssct,
            file = "P4_CellTypeASW_meanct.txt", #change
            quote =F,
            sep = "\t")
write.table(Medssct,
            file = "P4_CellTypeASW_medct.txt", #change
            quote =F,
            sep = "\t")


#----
#== CNB davies_bouldin_score ==#
setwd("/public/home/wutong/Project/PJ3_HumanLiver_scRNAseq_Databse/Result/Figure2/F2_last/F2_davies_bouldin_score/")
use_condaenv("/public/home/wutong/Software/miniconda3/envs/R4.1.1/bin/python")
py_module_available("sklearn")
repl_python()
import sklearn
from sklearn import metrics
from sklearn.metrics import davies_bouldin_score
SSlist_All = list()
SSlist_Center = list()
SSlist_Normal = list()
SSlist_Border = list()
for n in range(15):
  tmpSS_All = davies_bouldin_score(r.DataPCAlist_All, r.DataMetalist_All[:,n])
  tmpSS_Center = davies_bouldin_score(r.DataPCAlist_Center, r.DataMetalist_Center[:,n])
  tmpSS_Normal = davies_bouldin_score(r.DataPCAlist_Normal, r.DataMetalist_Normal[:,n])
  tmpSS_Border = davies_bouldin_score(r.DataPCAlist_Border, r.DataMetalist_Border[:,n])
  SSlist_All.append(tmpSS_All)
  SSlist_Center.append(tmpSS_Center)
  SSlist_Normal.append(tmpSS_Normal)
  SSlist_Border.append(tmpSS_Border)
quit

ss_list <- list()
ss_list$All <- data.frame(py$SSlist_All)
ss_list$Center <- data.frame(py$SSlist_Center) #change
ss_list$Normal <- data.frame(py$SSlist_Normal)
ss_list$Border <- data.frame(py$SSlist_Border)

sslist2 <- list()
meansslist <- list()
medsslist <- list()
#meanssctlist <- list()
#medssctlist <- list()
for(z in 1:length(ss_list)){
  name <- head(scRNAlist[[z]]@meta.data$CellClass,1)
  tmp_ss <- ss_list[[z]] %>%
    set_names(str_c("res", seq(0.1,1.5,0.1), sep = "")) %>%
    t() %>%
    data.frame() %>%
    mutate(Dataset = rep(name)) 
  tmp_ss <- tmp_ss %>% mutate(Resolution = row.names(tmp_ss))
  rownames(tmp_ss) <- NULL
  colnames(tmp_ss) <- c("davies_bouldin_score","CellType","Resolution") 
  mean_ss <- mean(tmp_ss$davies_bouldin_score) 
  med_ss <- median(tmp_ss$davies_bouldin_score)
  sslist2[[z]] <- tmp_ss
  meansslist[[z]] <- mean_ss
  medsslist[[z]] <- med_ss
}
Allss <- do.call(rbind,sslist2)
Meanss <- do.call(rbind,meansslist)
Medss <- do.call(rbind,medsslist)


write.table(Allss,
            file = "P2_davies_bouldin_score.txt", #change
            quote =F,
            sep = "\t")
write.table(Meanss,
            file = "P2_davies_bouldin_score_mean.txt", #change
            quote =F,
            sep = "\t")
write.table(Medss,
            file = "P2_davies_bouldin_score_med.txt", #change
            quote =F,
            sep = "\t")

#Plot ASW & DBI
Allss$CellType <- as.factor(Allss$CellType)
Allss$CellClass <- factor(Allss$CellType, levels = c("All","Core","Middle","Marginal"))
P_F3A <- ggplot(Allss, aes(CellClass, ASW))+
  geom_boxplot(aes(fill = CellClass, color = CellClass), 
               alpha = 0.5, 
               outlier.shape = NA,
               width = 0.5,
               position = position_dodge(width = 0.6)) +
  #geom_jitter(aes(SampleName, JACCARD, color = SampleQuality), width = 0.05) +
  geom_dotplot(aes(CellClass, ASW, fill=CellClass, color = CellClass),
               binaxis = 'y',
               stackdir = 'center',
               dotsize = 0.4,
               position = position_dodge(width = 0.6)) +
  # scale_fill_manual(values = c("#3675b6","#86c4e7", "#bcbdbf"))+
  # scale_color_manual(values = c("#3675b6","#86c4e7", "#bcbdbf"))+
  scale_fill_manual(values = c("#F3C861","#e98307","#91BF6E","#b8c8db")) +
  scale_color_manual(values = c("#F3C861","#e98307","#91BF6E","#b8c8db")) +
  theme_bw() +
  labs(y = "ASW", x = "Cell Class") +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 10,colour = "black"),
        #legend.position = "none",
        axis.title.x = element_text(size = 10,colour = "black"), #X轴标题大小
        axis.title.y = element_text(size = 10,colour = "black"), #Y轴标题大小
        axis.text.x = element_text(size = 10,colour = "black"), #X轴刻度标签大小
        axis.text.y = element_text(size = 10,colour = "black")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(panel.grid = element_blank()) +
  ylim(0.05, 0.32) 
#  ylim(0,1)

stat.test <-  Allss %>%
  t_test(ASW~CellClass, paired = T) 

stat.test2 <- stat.test %>% filter(group1 == "Core" | group2 == "Core")
stat.test2 <- stat.test2 %>% add_xy_position() %>%
  mutate(xmin = c(1, 2.05, 2.05)) %>%
  mutate(xmax = c(1.95, 3, 4)) %>%
  mutate(y.position = c(0.315,0.315,0.32))#添加XY位置
P_F3A <- P_F3A + stat_pvalue_manual(stat.test2, label = "p.adj.signif",
                                    bracket.size = 0.5, # 粗细
                                    tip.length = 0.01,
                                    label.size = 3)
ggsave("P2_ASW.pdf", 
       plot = P_F3A, width = 3.8,height = 2.7)
saveRDS(stat.test, "P2_ASW_test.rds")

