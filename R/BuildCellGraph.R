#' Build cell–cell distance matrix from PCA embeddings
#' This function extracts PCA embeddings from a Seurat object 
#' and computes the Euclidean distance between all cells to 
#' generate a distance matrix.
#'
#' @param seurat_obj A Seurat object containing dimensionality reduction results.
#' @param reduction A string specifying which reduction to use (default is "pca").
#'
#' @return A cell-cell distance matrix.
#' @export
#'
#' @importFrom Seurat Embeddings
#' @importFrom stats dist
#'
#' @examples
#' \dontrun{
#' dists <- BuildCellGraph(seurat_obj, reduction = "pca")
#' }

BuildCellGraph <- function(seurat_obj, reduction = "pca") {
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object.")
  }
  if (!(reduction %in% names(seurat_obj@reductions))) {
    stop(paste("Reduction", reduction, "not found in seurat_obj."))
  }
  pca_mat <- Seurat::Embeddings(seurat_obj, reduction = reduction)
  d <- as.matrix(dist(pca_mat))
  return(d)
}
