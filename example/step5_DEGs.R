library(pacman)
pacman::p_load(ggplot2,dplyr,tidyr,Matrix,RobustRankAggreg,
               stringr,Seurat,cowplot, clustree,
               tibble,patchwork,SCpubr,ROGUE,
               Scillus,magrittr,tidyverse,reticulate,
               ggsci,reshape2,scDblFinder,SingleR,celldex)
setwd("/public/home/wutong/test") #change
use_condaenv("/public/home/wutong/Software/miniconda3/envs/R4.1.1/bin/python")
py_module_available("sklearn") 

#----
#== Find all markers ==#
# Load single-cell data (all cells)
bf_scRNA <- readRDS(rds1)
# Load single-cell data (core cells only)
aft_scRNA <- readRDS(rds2)

bf_marker_wilcox <- list()
bf_marker_bimod <- list()
bf_marker_roc <- list()
bf_marker_t <- list()
bf_marker_negbinom <- list()
bf_marker_poisson <- list()
bf_marker_LR <- list()
bf_marker_MAST <- list()
bf_marker_DESeq2 <- list()

aft_marker_wilcox <- list()
aft_marker_bimod <- list()
aft_marker_roc <- list()
aft_marker_t <- list()
aft_marker_negbinom <- list()
aft_marker_poisson <- list()
aft_marker_LR <- list()
aft_marker_MAST <- list()
aft_marker_DESeq2 <- list()

for (i in 1:4) {
  bf_marker_wilcox[[i]] <- FindAllMarkers(bf_scRNA[[i]], assay = "RNA",test.use = "wilcox")
  bf_marker_bimod[[i]] <- FindAllMarkers(bf_scRNA[[i]], assay = "RNA",test.use = "bimod")
  bf_marker_roc[[i]] <- FindAllMarkers(bf_scRNA[[i]], assay = "RNA",test.use = "roc")
  bf_marker_t[[i]] <- FindAllMarkers(bf_scRNA[[i]], assay = "RNA",test.use = "t")
  bf_marker_negbinom[[i]] <- FindAllMarkers(bf_scRNA[[i]], assay = "RNA",test.use = "negbinom")
  bf_marker_poisson[[i]] <- FindAllMarkers(bf_scRNA[[i]], assay = "RNA",test.use = "poisson")
  bf_marker_LR[[i]] <- FindAllMarkers(bf_scRNA[[i]], assay = "RNA",test.use = "LR")
  bf_marker_MAST[[i]] <- FindAllMarkers(bf_scRNA[[i]], assay = "RNA",test.use = "MAST")
  bf_marker_DESeq2[[i]] <- FindAllMarkers(bf_scRNA[[i]], assay = "RNA",test.use = "DESeq2")
  #
  aft_marker_wilcox[[i]] <- FindAllMarkers(aft_scRNA[[i]], assay = "RNA",test.use = "wilcox")
  aft_marker_bimod[[i]] <- FindAllMarkers(aft_scRNA[[i]], assay = "RNA",test.use = "bimod")
  aft_marker_roc[[i]] <- FindAllMarkers(aft_scRNA[[i]], assay = "RNA",test.use = "roc")
  aft_marker_t[[i]] <- FindAllMarkers(aft_scRNA[[i]], assay = "RNA",test.use = "t")
  aft_marker_negbinom[[i]] <- FindAllMarkers(aft_scRNA[[i]], assay = "RNA",test.use = "negbinom")
  aft_marker_poisson[[i]] <- FindAllMarkers(aft_scRNA[[i]], assay = "RNA",test.use = "poisson")
  aft_marker_LR[[i]] <- FindAllMarkers(aft_scRNA[[i]], assay = "RNA",test.use = "LR")
  aft_marker_MAST[[i]] <- FindAllMarkers(aft_scRNA[[i]], assay = "RNA",test.use = "MAST")
  aft_marker_DESeq2[[i]] <- FindAllMarkers(aft_scRNA[[i]], assay = "RNA",test.use = "DESeq2")
}

saveRDS(bf_marker_wilcox, "/public/home/wutong/Project/PJ3_HumanLiver_scRNAseq_Databse/Test/Eva_marker/bf/bf_marker_wilcox.rds")
saveRDS(bf_marker_bimod, "/public/home/wutong/Project/PJ3_HumanLiver_scRNAseq_Databse/Test/Eva_marker/bf/bf_marker_bimod.rds")
saveRDS(bf_marker_roc, "/public/home/wutong/Project/PJ3_HumanLiver_scRNAseq_Databse/Test/Eva_marker/bf/bf_marker_roc.rds")
saveRDS(bf_marker_t, "/public/home/wutong/Project/PJ3_HumanLiver_scRNAseq_Databse/Test/Eva_marker/bf/bf_marker_t.rds")
saveRDS(bf_marker_negbinom, "/public/home/wutong/Project/PJ3_HumanLiver_scRNAseq_Databse/Test/Eva_marker/bf/bf_marker_negbinom.rds")
saveRDS(bf_marker_poisson, "/public/home/wutong/Project/PJ3_HumanLiver_scRNAseq_Databse/Test/Eva_marker/bf/bf_marker_poisson.rds")
saveRDS(bf_marker_LR, "/public/home/wutong/Project/PJ3_HumanLiver_scRNAseq_Databse/Test/Eva_marker/bf/bf_marker_LR.rds")
saveRDS(bf_marker_MAST, "/public/home/wutong/Project/PJ3_HumanLiver_scRNAseq_Databse/Test/Eva_marker/bf/bf_marker_MAST.rds")
saveRDS(bf_marker_DESeq2, "/public/home/wutong/Project/PJ3_HumanLiver_scRNAseq_Databse/Test/Eva_marker/bf/bf_marker_DESeq2.rds")

saveRDS(aft_marker_wilcox, "/public/home/wutong/Project/PJ3_HumanLiver_scRNAseq_Databse/Test/Eva_marker/bf/aft_marker_wilcox.rds")
saveRDS(aft_marker_bimod, "/public/home/wutong/Project/PJ3_HumanLiver_scRNAseq_Databse/Test/Eva_marker/bf/aft_marker_bimod.rds")
saveRDS(aft_marker_roc, "/public/home/wutong/Project/PJ3_HumanLiver_scRNAseq_Databse/Test/Eva_marker/bf/aft_marker_roc.rds")
saveRDS(aft_marker_t, "/public/home/wutong/Project/PJ3_HumanLiver_scRNAseq_Databse/Test/Eva_marker/bf/aft_marker_t.rds")
saveRDS(aft_marker_negbinom, "/public/home/wutong/Project/PJ3_HumanLiver_scRNAseq_Databse/Test/Eva_marker/bf/aft_marker_negbinom.rds")
saveRDS(aft_marker_poisson, "/public/home/wutong/Project/PJ3_HumanLiver_scRNAseq_Databse/Test/Eva_marker/bf/aft_marker_poisson.rds")
saveRDS(aft_marker_LR, "/public/home/wutong/Project/PJ3_HumanLiver_scRNAseq_Databse/Test/Eva_marker/bf/aft_marker_LR.rds")
saveRDS(aft_marker_MAST, "/public/home/wutong/Project/PJ3_HumanLiver_scRNAseq_Databse/Test/Eva_marker/bf/aft_marker_MAST.rds")
saveRDS(aft_marker_DESeq2, "/public/home/wutong/Project/PJ3_HumanLiver_scRNAseq_Databse/Test/Eva_marker/bf/aft_marker_DESeq2.rds")

#----
#== DEG rank ==#
files<-dir(path = "/public/home/wutong/Project/PJ3_HumanLiver_scRNAseq_Databse/Result/step6_DE/V5/R1_Dataset_DE/DE_result/Center",
           full.names = T,
           pattern = ".rds")
marker1 <- map(files,readRDS)
markerW <- readRDS("/public/home/wutong/Project/PJ3_HumanLiver_scRNAseq_Databse/Result/step6_DE/V5/R1_Dataset_DE/DE_result/Center/marker_wilcox.rds")
marker2 <- list()
marker3 <- list()
marker4 <- list()
marker5 <- list()
marker6 <- list()
#df1<-reduce(df,inner_join)
#filter DEGs
for ( i in 1:10) {
  for ( n in 1:length(marker1)) {
    #基因-细胞簇列表，不同检验方法中取交集
    marker2[[n]] <- marker1[[n]][[i]] %>% 
      mutate(gene_cluster = str_c(gene, cluster, sep = ":")) %>%
      select(gene_cluster) %>% 
      as.data.frame()
  }
  #single cell DEGs overlap in different test
  marker3[[i]] <- purrr::reduce(marker2, inner_join)
  marker4[[i]] <- markerW[[i]] %>%
    mutate(gene_cluster = str_c(gene, cluster, sep = ":")) %>%
    inner_join(marker3[[i]], by = "gene_cluster") %>%
    filter(avg_log2FC >= 0.5) %>%
    filter(p_val_adj <= 0.01) %>%
    filter(pct.1 >= 0.2) %>%
    mutate(pv2 = if_else(p_val_adj <= 1e-310, 1e-310, p_val_adj)) %>%
    mutate(pv3 = -log10(pv2)) %>%
    mutate(rank_score = avg_log2FC * pv3) %>%
    mutate(dataset_name = str_c("D", i, sep = ""))#单个数据集中所有细胞类型marker rank score list
  marker4[[i]]$rank_score <- as.numeric(marker4[[i]]$rank_score)
}
marker5 <- do.call(rbind, marker4)
marker6 <- marker5 %>%
  as_tibble() %>%
  group_split(cluster) 

rank_ag_list <- list()
rank_ag_list2 <- list()
for (m in 1:length(marker6)) {
  clustername = marker6[[m]]$cluster[1] %>% as.character()
  #marker6[[m]]$rank_score <- as.integer(marker6[[m]]$rank_score)
  tmp <- marker6[[m]] %>% as.data.frame()
  tmp <- tmp %>% group_split(dataset_name)
  tmp_rank_list <- list()
  for ( x in 1:length(tmp)) {
    tmp[[x]] <- tmp[[x]] %>% 
      as.data.frame() %>%
      #按照avg_log2FC * -log10(pv2)数值从大到小排序
      #log2FC越大，Pv越小，排名越靠前
      arrange(-rank_score) 
    tmp_rank_list[[x]] <- tmp[[x]]$gene
  }
  # Based on the product of vg_log2FC and -log10(pv2)
  rank_ag_list[[m]] = aggregateRanks(tmp_rank_list)
  freq = as.data.frame(table(unlist(tmp_rank_list))) 
  freq <- freq %>%
    dplyr::rename(Gene = Var1)
  # Add the frequency of each gene's occurrence across the 10 datasets to the aggregated ranking results
  rank_ag_list[[m]] <- rank_ag_list[[m]] %>% 
    dplyr::rename(Gene = Name) %>%
    mutate(cluster = rep(clustername)) %>% left_join(freq, by = "Gene")
  tmprank <- rank_ag_list[[m]]
  tmprank <- tmprank %>% group_split(Freq)
  tmprank2 <- list()
  for ( x in 1:length(tmprank)) {
    tmprank2[[x]] <- tmprank[[x]] %>%
      # Compute a reverse order index for each group
      mutate(rank1_freq = rep(length(tmprank)+1-x)) %>% 
      arrange(Score) %>%
      mutate(rank2_score = 1:nrow(tmprank[[x]]))
  }
  rank_ag_list2[[m]] <- do.call(rbind, tmprank2)
  rank_ag_list2[[m]] <- rank_ag_list2[[m]] %>% 
    mutate(rank_score = rank1_freq + (rank2_score / 10000)) %>%
    arrange(rank_score) %>%
    mutate(Last_rank = 1:nrow(rank_ag_list2[[m]]))
}
saveRDS(rank_ag_list2, "/public/home/wutong/Project/PJ3_HumanLiver_scRNAseq_Databse/Result/step6_DE/V5/R1_Dataset_DE/Eva_result/Center_DE_rank_PctMinus0.rds")

#----
#== Marker ROC ==#
#----
library(ROCR)
library(ggplot2)
library(dplyr)

#读入文件
# Load a list of re-integrated RDS objects from four populations (all cells included)
scRNAlist <- readRDS("All_core_inter_marg.rds")
# Load 145 gold-standard markers collected from the literature
marker <- readRDS("AnnoMarker145.rds")
setwd("/your_path/")

# Define the cell types of interest and their corresponding signature genes
cellname <- "lrNK" #
gene <- "XCL1"
NM <- str_c(cellname,gene,"ROC_Curve.pdf",sep = "_")

roc_data_list <- list() # Container for ROC data of all genes

for (i in 1:length(scRNAlist)) {
  RNA <- scRNAlist[[i]]@assays$RNA@counts
  celltype <- scRNAlist[[i]]@meta.data %>%
    select(CellType,Cell1)
  celltype2 <- scRNAlist[[i]]@meta.data %>% select(CellType) %>% distinct()
  # Retrieve gene expression data
  exp <- RNA[gene,]
  exp2 <- as.data.frame(exp) 
  exp2 <- exp2 %>% mutate(Cell1 = rownames(exp2))
  exp3 <- merge(exp2, celltype, "Cell1")
  exp3 <- exp3 %>% mutate(CellType2 = CellType)
  exp3$CellType2 <- ifelse(is.na(exp3$CellType), 1, ifelse(exp3$CellType == cellname, 0, 1))
  colnames(exp3) <- c("Cell1", "exp","CellType", "CellType2")
  # AUC
  prediction_gene <- NULL
  prediction_gene <- prediction(predictions = exp3$exp, labels = exp3$CellType2, label.ordering = 1:0)
  perf_gene <- performance(prediction.obj = prediction_gene, measure = "tpr", x.measure = "fpr")
  auc_gene <- performance(prediction.obj = prediction_gene, measure = "auc")
  auc_value <- round(auc_gene@y.values[[1]], digits = 3)
  
  # ROC
  fpr_gene <- perf_gene@x.values[[1]]
  tpr_gene <- perf_gene@y.values[[1]]
  
  # Create ROC dataframe
  roc_data_gene <- data.frame(FPR = fpr_gene, TPR = tpr_gene, Gene = gene, AUC = auc_value)
  
  #添加Cell Class类型
  cellclass <- scRNAlist[[i]]$CellClass[[1]]
  roc_data_gene <- roc_data_gene %>% mutate(CellClass = cellclass)
  roc_data_list[[i]] <- roc_data_gene
}
roc_data <- do.call(rbind, roc_data_list)
roc_data$CellClass <- gsub("Center","Core",roc_data$CellClass)
roc_data$CellClass <- gsub("Normal","Middle",roc_data$CellClass)
roc_data$CellClass <- gsub("Border","Marginal",roc_data$CellClass)
saveRDS(roc_data, str_c(cellname,gene,"ROC_Curve.rds",sep = "_"))

# Annotate cells with their class type
All_AUC <- roc_data %>% filter(CellClass == "All") %>% select(AUC) %>% head(1)
Core_AUC <- roc_data %>% filter(CellClass == "Core") %>% select(AUC) %>% head(1)
Middle_AUC <- roc_data %>% filter(CellClass == "Middle") %>% select(AUC) %>% head(1)
Marginal_AUC <- roc_data %>% filter(CellClass == "Marginal") %>% select(AUC) %>% head(1)

colors <- c("All" = "#F3C861", 
            "Core" = "#e98307", 
            "Middle" = "#91BF6E", 
            "Marginal" = "#b8c8db")

# Plot ROC curves for the selected markers
test <- ggplot(roc_data, aes(x = FPR, y = TPR, color = CellClass, group = CellClass)) +
  geom_line(size = 0.8) +  # 根据 CellClass 绘制不同颜色的线条
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#c4c4c4", size = 0.8) +  # 添加对角线
  labs(x = "False positive rate", y = "True positive rate") +
  scale_color_manual(values = colors) +  # 为不同 CellClass 指定颜色
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        axis.text.y = element_text(angle = 0, color = "black", size = 14),
        axis.text.x = element_text(angle = 0, color = "black", size = 14),
        axis.title.y = element_text(size = 14),  
        axis.title.x = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 14),
        legend.position = "right") +
  # Add annotation based on AUC values
  annotate("text", x = 0.4, y = 0.24, label = paste("All AUC = ",All_AUC), color = "#F3C861", size = 4, hjust = 0) +
  annotate("text", x = 0.4, y = 0.17, label = paste("Core AUC = ",Core_AUC), color = "#e98307", size = 4, hjust = 0) +
  annotate("text", x = 0.4, y = 0.1, label = paste("Middle AUC = ",Middle_AUC), color = "#91BF6E", size = 4, hjust = 0) +
  annotate("text", x = 0.4, y = 0.03, label = paste("Marginal AUC = ",Marginal_AUC), color = "#b8c8db", size = 4, hjust = 0)

# Save the generated image
ggsave(NM, plot = test, width = 4.6, height = 3.3, limitsize = FALSE, dpi = 6000)

#----
#== AUROC ==#
#gene expression matrix
scRNA <- readRDS("/public/home/wutong/Project/PJ3_HumanLiver_scRNAseq_Databse/Result/step8_celltypist_classifer/ValidateDatasets/D2_GSE174748/h-livernew_celltype_V2.rds")
RNA <- scRNA@assays$RNA@counts
celltype <- scRNA@meta.data %>% 
  select(CellType,Cell1)  %>%
  dplyr::rename(CellType = celltype)
celltype2 <- scRNA@meta.data %>% select(celltype) %>% distinct()

# filter row sum = 0
library(Matrix)
row_sums <- rowSums(RNA)
RNA_filter <- RNA[row_sums != 0, ]

# gene list
## Paper marker gene
CellTaxonomy_marker <- readRDS("/public/home/wutong/Project/PJ3_HumanLiver_scRNAseq_Databse/Prepare/Database_Marker/CellTaxonomy_marker/Cell_Taxonomy_HumanLiver.rds")
#CellTaxonomy_marker <- CellTaxonomy_marker %>% filter(Tissue_standard != "Blood")
test1 <- CellTaxonomy_marker %>% 
  filter(Cluster %in% c("B cell","CD8+ T","Plasma cell","Mast cell","pDC")) %>% 
  filter(!is.na(Condition))
test2 <- CellTaxonomy_marker %>% 
  filter(Cluster == "pDC") %>% 
  filter(Tissue_standard != "Blood")
CellTaxonomy_marker_tmp <- CellTaxonomy_marker %>% 
  filter(!(Cluster %in% c("B cell","CD8+ T","Plasma cell","Mast cell","pDC")))
CellTaxonomy_marker <- rbind(CellTaxonomy_marker_tmp,test1,test2)
CellTaxonomy_marker2 <- CellTaxonomy_marker %>% 
  select(Cluster, Cell_Marker) %>% 
  mutate(cluster_gene = str_c(Cluster, Cell_Marker, sep = "_")) %>%
  dplyr::rename(tmpname = Cluster, tmpname2 = Cell_Marker) 

papermarker <- distinct(CellTaxonomy_marker2) # 去重
colnames(papermarker) <- c("cluster","gene","cluster_gene")
gene1 <- papermarker$tmpname2
gene2 <- rownames(RNA)
gene3 <- intersect(gene1,gene2) # Intersect with genes in the RNA matrix to prevent errors
papermarker <- papermarker %>% 
  filter(gene %in% gene3) %>% 
  filter(cluster %in% celltype2$celltype)

Genelist <- as.list(papermarker)
Genelist2 <- papermarker %>% group_split(cluster)

## annotation marker gene 
Gene <- readRDS("/public/home/wutong/Project/PJ3_HumanLiver_scRNAseq_Databse/Result/Figure2/F2_last/F2_F_AUC/F2_F_AnnoMarker.rds")
colnames(Gene) <- c("cluster","gene","cluster_gene")
gene1 <- Gene$gene 
gene2 <- rownames(RNA)
gene3 <- intersect(gene1,gene2) 
Gene <- Gene %>% filter(gene %in% gene3) %>% 
  filter(cluster %in% celltype2$celltype)
Genelist <- as.list(Gene)
Genelist2 <- Gene %>% group_split(cluster)

## Rank DE gene 
Center_DE <- readRDS("/public/home/wutong/Project/PJ3_HumanLiver_scRNAseq_Databse/Result/step6_DE/V5/R1_Dataset_DE/Eva_result/Center_DE_rank_PctMinus0.rds")
DEtop <- do.call(rbind, Center_DE)

gene1 <- DEtop$Gene
gene3 <- intersect(gene1,gene2)
DEtop <- DEtop %>% filter(Gene %in% gene3)

DErank_list <- list()
for (i in 1:length(Genelist2)) {
  tmp <- Genelist2[[i]]
  CellType <- tmp$cluster[1]
  num <- dim(tmp)[1][[1]]
  DErank_list[[i]] <- DEtop %>% filter(cluster == CellType) %>% slice_min(order_by = Last_rank, n = num)
}
DErank <- do.call(rbind, DErank_list) %>% 
  select(cluster, Gene) %>% 
  dplyr::rename(gene = Gene)
Genelist <- as.list(DErank)

#== Gene AUC ==#
group = "All"
AUClist <- list()
for (i in 1:length(Genelist[[1]])) {
  cellname <- Genelist[[1]][[i]]
  genename <- Genelist[[2]][[i]]
  print(genename)
  exp <- RNA[genename,]
  exp2 <- as.data.frame(exp) 
  exp2 <- exp2 %>% mutate(Cell1 = rownames(exp2))
  exp3 <- merge(exp2, celltype, "Cell1")
  exp3 <- exp3 %>% mutate(CellType2 = CellType)
  exp3$CellType2 <- ifelse(exp3$CellType == cellname, 0, 1)
  # AUC: set label.ordering = 0:1 (0 = positive class, 1 = negative class)
  prediction.use <- prediction(predictions = exp3$exp, labels = exp3$CellType2, label.ordering = 1:0)
  perf.use <- performance(prediction.obj = prediction.use, measure = "auc")
  auc.use <- round(x = perf.use@y.values[[1]], digits = 3)
  AUClist[[i]] <- cbind(genename, auc.use, cellname, group)
}
AUC <- do.call(rbind, AUClist)
AnnotateMarker_AUC_All <- as.data.frame(AUC)

#----
#== Plot ==#
AllAUC <- rbind(RankDE_AUC, DatabaseMarker_AUC)
colnames(AllAUC) <- c("genename","myAUC","cluster","CellClass")
AllAUC$CellClass <- factor(AllAUC$CellClass, levels =c("AnnotateMarker","RankDE"))
AllAUC$myAUC <- as.numeric(AllAUC$myAUC) 

#cluster order
RankDE_AUC$auc.use <- as.numeric(RankDE_AUC$auc.use)
orderlist <- RankDE_AUC %>% group_split(cellname) 

orderlist2 <- list()
for (i in 1:length(orderlist)) {
  orderlist2[[i]] <- cbind(as.character(orderlist[[i]]$cellname[[1]]), median(orderlist[[i]]$auc.use))
}
orderlist3 <- do.call(rbind, orderlist2) %>% as.data.frame()
orderlist3$V2 <- as.numeric(orderlist3$V2)
orderlist4 <- orderlist3 %>% arrange(-V2) %>% select(V1) 

AllAUC$cluster <- factor(AllAUC$cluster, levels = orderlist4$V1)

F4_C <- ggplot(AllAUC, aes(cluster, myAUC))+
  geom_boxplot(aes(fill = CellClass, color = CellClass), 
               alpha = 0.5, 
               outlier.shape = NA,
               width = 0.5,
               position = position_dodge(width = 0.6)) +
  #geom_jitter(aes(SampleName, JACCARD, color = SampleQuality), width = 0.05) +
  # geom_dotplot(aes(cluster, myAUC, fill=CellClass,color = CellClass),
  #              binaxis = 'y',
  #              stackdir = 'center',
  #              dotsize = 0.2,
  #              position = position_dodge(width = 0.6)) +
  scale_fill_manual(values = c("#eac590","#225f9b"))+
  scale_color_manual(values = c("#eac590","#225f9b"))+
  scale_fill_manual(values = c("#eac590","#225f9b")) +
  scale_color_manual(values = c("#eac590","#225f9b")) +
  theme_bw() +
  labs(y = "AUC", x = "CellType") +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 10),
        #legend.position = "none",
        axis.title.x = element_text(size = 10,colour = "black"), 
        axis.title.y = element_text(size = 10,colour = "black"), 
        axis.text.x = element_text(size = 10,colour = "black"), 
        axis.text.y = element_text(size = 10,colour = "black")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(panel.grid = element_blank()) +
  ylim(0.4, 1.03)


stat.test <- AllAUC %>% 
  #data_new1[which(data_new1$group_number == 2),] %>%
  group_by(cluster) %>%
  pairwise_t_test(
    myAUC ~ CellClass, 
    #paired = FALSE, 
    p.adjust.method = "fdr") %>%
  add_xy_position(x = "cluster")


stat2 <- stat.test %>% 
  mutate(y.position = y.position - 0.02) 


F4_C <- F4_C + stat_pvalue_manual(stat.test, label = "p.adj.signif",
                                  bracket.size = 0.5, # 粗细
                                  tip.length = 0.01,
                                  label.size = 2.8)


ggsave("F4_C_AnnoMarker_RankDE_AUC.pdf", 
       plot = F4_C, width = 12,height = 5)
saveRDS(stat.test, "F4_C_AnnoMarker_AUC_Ttest.rds")
saveRDS(AllAUC, "F4_C_AnnoMarker_AUC.rds")
saveRDS(DErank, "F4_C_AnnoMarker_DErank.rds")
saveRDS(Gene, "F4_C_AnnoMarker_.rds")


