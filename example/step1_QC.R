#----
#== Packages ==#
library(pacman)
pacman::p_load(ggplot2,dplyr,tidyr,Matrix,Rmisc,
               stringr,Seurat,cowplot,clustree,
               tibble,patchwork,SCpubr,ROGUE,scuttle,
               Scillus,magrittr,tidyverse,reticulate,
               ggsci,reshape2,scDblFinder,SingleR,celldex)
setwd("/public/home/wutong/test") #change
use_condaenv("/public/home/wutong/Software/miniconda3/envs/R4.1.1/bin/python")
py_module_available("sklearn") #check

#----
#== Part2 QC ==#
# Load a list of single-sample RDS objects (each element corresponds to one sample)
scRNAlist <- readRDS("sample_list.rds")
DbState <- list()
nFeatureMed <- list()
nCountMed <- list()
QC_state <- list()
QC_FeatureScatter <- list()
QCnF <- list()
QCnC <- list()
s.genes=Seurat::cc.genes.updated.2019$s.genes
g2m.genes=Seurat::cc.genes.updated.2019$g2m.genes

for ( i in 1:length(tmp_scRNAlist)) {
  tmp_scRNAlist[[i]] <- tmp_scRNAlist[[i]][!grepl("^MT-", rownames(tmp_scRNAlist[[i]])), ]
  tmp_scRNAlist[[i]] <- tmp_scRNAlist[[i]][!grepl("^RP[SL]", rownames(tmp_scRNAlist[[i]])), ]
  tmp_scRNAlist[[i]] <- tmp_scRNAlist[[i]][!grepl("^ERCC", rownames(tmp_scRNAlist[[i]])), ]
  testnf <- isOutlier(tmp_scRNAlist[[i]]$nFeature_RNA, 
                      nmads = 3, 
                      log=TRUE,
                      type="both")
  testnc <- isOutlier(tmp_scRNAlist[[i]]$nCount_RNA, 
                      nmads = 3, 
                      log=TRUE,
                      type="both")
  lowernf <- attr(testnf, "thresholds")[1] %>% as.numeric()
  highernf <- attr(testnf, "thresholds")[2] %>% as.numeric()
  lowernc <- attr(testnc, "thresholds")[1] %>% as.numeric()
  highernc <- attr(testnc, "thresholds")[2] %>% as.numeric()
  QCnF[[i]] <- attr(testnf, "thresholds")
  QCnC[[i]] <- attr(testnc, "thresholds")
  tmp_scRNAlist[[i]] <- subset(tmp_scRNAlist[[i]], subset = nFeature_RNA > lowernf & nFeature_RNA < highernf & nCount_RNA > lowernc &  nCount_RNA < highernc & percent_mito < 20)
  scRNAlist[[i]] <- tmp_scRNAlist[[i]]
  tmp_scRNAlist[[i]] <- scDblFinder(as.SingleCellExperiment(tmp_scRNAlist[[i]]))
  scRNAlist[[i]]$scDblFinder.score <- tmp_scRNAlist[[i]]$scDblFinder.score
  scRNAlist[[i]]$scDblFinder.class <- tmp_scRNAlist[[i]]$scDblFinder.class
  nCount <- scRNAlist[[i]]@meta.data$nCount_RNA
  nCountMed[[i]] <- median(nCount)
  nFeature <- scRNAlist[[i]]@meta.data$nFeature_RNA
  nFeatureMed[[i]] <- median(nFeature)
  DbState[[i]] <- table(tmp_scRNAlist[[i]]$scDblFinder.class)
  scRNAlist[[i]] <- subset(scRNAlist[[i]],
                           subset = scDblFinder.class=='singlet')
  scRNAlist[[i]] <- NormalizeData(scRNAlist[[i]])
  scRNAlist[[i]] <- FindVariableFeatures(scRNAlist[[i]], selection.method = "vst")
  scRNAlist[[i]] <- CellCycleScoring(scRNAlist[[i]],
                                     s.features = s.genes,
                                     g2m.features = g2m.genes,
                                     set.ident = TRUE)
}

scRNA_merger_aftQC <- merge(scRNAlist[[1]],
                           list(scRNAlist[[2]],
                                scRNAlist[[3]],
                                scRNAlist[[4]],
                                scRNAlist[[5]])) 
#QCtable
Singlet <- list()
Doublet <- list()
DbRate <- list()
LQCnF <- list()
HQCnF <- list()
LQCnC <- list()
HQCnC <- list()
for(i in 1:(length(QCnF))){LQCnF[[i]] <- QCnF[[i]][["lower"]]}
for(i in 1:(length(QCnF))){HQCnF[[i]] <- QCnF[[i]][["higher"]]}
for(i in 1:(length(QCnC))){LQCnC[[i]] <- QCnC[[i]][["lower"]]}
for(i in 1:(length(QCnC))){HQCnC[[i]] <- QCnC[[i]][["higher"]]}
for(i in 1:(length(DbState))){Singlet[[i]] <- DbState[[i]][["singlet"]]}
for(i in 1:(length(DbState))){Doublet[[i]] <- DbState[[i]][["doublet"]]}
for(i in 1:(length(DbState))){DbRate[[i]] <- DbState[[i]][["doublet"]]/(DbState[[i]][["doublet"]]+DbState[[i]][["singlet"]])}
Ncell_bf <- table(scRNA_merger_bfQC$batch) %>%
  as.data.frame() %>%
  tidygraph::rename(SampleName = Var1, Ncell_bfQC = Freq)
Ncell_aft <- table(scRNA_merger_aftQC$batch) %>%
  as.data.frame() %>%
  select(-Var1) %>%
  tidygraph::rename(Ncell_aftQC = Freq)
state_table <- cbind(nFeatureMed,nCountMed,LQCnF,HQCnF,LQCnC,HQCnC,Singlet,Doublet,DbRate)
state_table <-  state_table %>% 
  as.data.frame() %>%
  mutate(SN = Ncell_bf$SampleName, Ncell_bfQC = Ncell_bf$Ncell_bfQC, Ncell_aftQC = Ncell_aft$Ncell_aftQC) 
state_table <- data.frame(lapply(state_table, as.character), stringsAsFactors=FALSE)

#Plot
aftqc_phisnc <- list()
aftqc_phisnf <- list()
for ( i in 1:length(tmp_scRNAlist)) {
  aftqc_phisnf[[i]] <- ggplot(data_frame(scRNAlist[[i]]@meta.data), aes(nFeature_RNA)) + geom_histogram(binwidth = 100, bins = 300) + ggtitle(str_c("s",i,"_nFeature"))
  aftqc_phisnc[[i]] <- ggplot(data_frame(scRNAlist[[i]]@meta.data), aes(nCount_RNA)) + geom_histogram(binwidth = 100, bins = 300) + ggtitle(str_c("s",i,"_nCount"))
}
aftqc_p_his_nf <- wrap_plots(aftqc_phisnf, byrow = T, ncol = 3) #根据样本数量调整列数
aftqc_p_his_nc <- wrap_plots(aftqc_phisnc, byrow = T, ncol = 3) #根据样本数量调整列数

aftqc_p_vil <- VlnPlot(scRNA_merger_aftQC, 
                 pt.size = 0.1,
                 group.by = "batch",
                 features = c("nFeature_RNA", "nCount_RNA", "percent_mito","percent_ribo","percent_hb"), 
                 ncol = 5)

saveRDS(scRNAlist, file = "/public/home/wutong/test/GSE115469/Pt2_aftQC_scRNAlist.rds") #改路径
saveRDS(scRNA_merger_aftQC, file = "/public/home/wutong/test/GSE115469/Pt2_aftQC_scRNA_merge.rds") 
write.table(state_table,"/public/home/wutong/test/GSE115469/Pt2_state_qc_table.xls", quote = F, sep = "\t", row.names = F)
ggsave("/public/home/wutong/test/GSE115469/Pt2_aftQC_his_nF.pdf", plot = aftqc_p_his_nf, width = 6, height = 4) #根据样本数调整输出长宽
ggsave("/public/home/wutong/test/GSE115469/Pt2_aftQC_his_nC.pdf", plot = aftqc_p_his_nc, width = 6, height = 4) 
ggsave("/public/home/wutong/test/GSE115469/Pt2_aftQC_vil.pdf", plot = aftqc_p_vil, width = 10, height = 4) 

