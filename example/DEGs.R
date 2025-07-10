library(optparse)
option_list <- list(
  make_option(c("-r", "--rds1"), type = "character", help = "rds1"),
  make_option(c("-R", "--rds2"), type = "character", help = "rds2")
)

opt <- parse_args(OptionParser(option_list=option_list))

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
bf_scRNA <- readRDS(opt$rds1)
aft_scRNA <- readRDS(opt$rds2)

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
