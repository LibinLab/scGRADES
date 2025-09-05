#' Perform clustering at multiple resolutions and return Seurat object
#' This function applies Seurat's graph-based clustering (FindClusters)
#' across a vector of resolution values and returns the Seurat object with
#' clustering results saved in its meta.data.
#' 
#' @title FindClusterAcrossRes
#' @description Perform clustering at multiple resolutions and return Seurat object
#' @name FindClusterAcrossRes
#' @param seurat_obj A Seurat object that has already undergone dimensionality reduction (e.g., PCA).
#' @param resolutions A numeric vector of resolution values. Default is seq(0.1, 1.5, 0.1).
#' @param assay The assay to use for clustering (e.g., "integrated" or "RNA"). Default is "integrated".
#'
#' @return A Seurat object with clustering results at all specified resolutions saved in meta.data.
#' @export
#'
#' @importFrom Seurat DefaultAssay<- FindNeighbors FindClusters 
#' @importFrom cli cli_process_start cli_process_done
#'
#' @examples
#' \dontrun{
#' seurat_obj <- FindClusterAcrossRes(seurat_obj, 
#'                                   resolutions = seq(0.2, 1.0, 0.2), 
#'                                   assay = "integrated")
#' }
   
FindClusterAcrossRes <- function(...) {
  id <- cli::cli_process_start("Building graph & clustering across resolutions ...")
  on.exit(cli::cli_process_done(id), add = TRUE)
  FindClusterAcrossRes_inter(...)
}                         

FindClusterAcrossRes_inter <- function(seurat_obj,
                                    resolutions = seq(0.1, 1.5, 0.1),
                                    assay = "integrated") {
  #require(Seurat)
  Seurat::DefaultAssay(seurat_obj) <- assay
  # 构建邻接图（默认用 assay_snn）
  seurat_obj <- Seurat::FindNeighbors(seurat_obj)
  for (res in resolutions) {
    seurat_obj <- Seurat::FindClusters(seurat_obj, resolution = res)
  }
  return(seurat_obj)
}
