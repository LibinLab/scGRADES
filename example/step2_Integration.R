library(Seurat)
library(ggplot2)
library(profmem)
library(rlog)
library(future)

# Load a list of RDS objects after QC for individual samples
dataset_merge.list=readRDS("dataset_merge_list_v4.rds")
features <- SelectIntegrationFeatures(object.list = dataset_merge.list, nfeatures = 5000)

rlog::log_info('FindIntegrationAnchors')

system.time({
  dataset_merge.anchors <- FindIntegrationAnchors(object.list = dataset_merge.list, anchor.features = features, reduction = "rpca")
  })

system.time({
  dataset_merge_integrated <- IntegrateData(anchorset = dataset_merge.anchors)
  })# output is seurat object

rlog::log_info('analysis of integrated object')
DefaultAssay(dataset_merge_integrated) <- "integrated"
dataset_merge_integrated <- ScaleData(dataset_merge_integrated)
dataset_merge_integrated <- RunPCA(dataset_merge_integrated, npcs = 30)
dataset_merge_integrated <- RunUMAP(dataset_merge_integrated, reduction = "pca", dims = 1:30)
dataset_merge_integrated <- FindNeighbors(dataset_merge_integrated, reduction = "pca", dims = 1:30)
dataset_merge_integrated <- FindClusters(dataset_merge_integrated, resolution = 1)

#save file
saveRDS(dataset_merge_integrated, file = "dataset_integrated_analyzed.rds") 

#visualisation
p1 <- DimPlot(dataset_merge_integrated, reduction = "umap", group.by = "samplebatch",raster=FALSE)
p2 <- DimPlot(dataset_merge_integrated, reduction = "umap", group.by = "sample",raster=FALSE)
p3 <- DimPlot(dataset_merge_integrated , reduction = "umap", label = TRUE, repel = TRUE,raster=FALSE)
ggsave("p1_test_res1.pdf", plot = p1+p2+p3, width = 45, height = 10)

