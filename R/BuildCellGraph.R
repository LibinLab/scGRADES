#' Build cell–cell distance matrix from PCA embeddings
#' This function extracts PCA embeddings from a Seurat object 
#' and computes the Euclidean distance between all cells to 
#' generate a distance matrix.
#'
#' @title BuildCellGraph
#' @description Build cell–cell distance matrix from PCA embeddings
#' @name BuildCellGraph
#' @param seurat_obj A Seurat object containing dimensionality reduction results.
#' @param reduction A string specifying which reduction to use (default is "pca").
#'
#' @return A cell-cell distance matrix.
#' @export
#'
#' @importFrom Seurat Embeddings 
#' @importFrom stats dist
#' @importFrom bigdist bigdist
#' @importFrom cli cli_process_start cli_process_done
#'
#' @examples
#' \dontrun{
#' dists <- BuildCellGraph(seurat_obj, reduction = "pca")
#' }

BuildCellGraph <- function(seurat_obj, reduction = "pca") {
  id <- cli::cli_process_start("Computing PCA-based distance matrix ...")
  on.exit(cli::cli_process_done(id), add = TRUE)
  BuildCellGraph_inter(seurat_obj, reduction = "pca")
}

BuildCellGraph_inter <- function(seurat_obj, reduction = "pca") {
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object.")
  }
  if (!(reduction %in% names(seurat_obj@reductions))) {
    stop(paste("Reduction", reduction, "not found in seurat_obj."))
  }
  pca_mat <- Seurat::Embeddings(seurat_obj, reduction = reduction)
  cells <- rownames(pca_mat)
  #d <- as.matrix(dist(pca_mat))
  fprefix <- tempfile(pattern = "cell_dist_", tmpdir = tempdir())
  d <- bigdist::bigdist(pca_mat, method = "euclidean", file = fprefix)
  return(d)
}

BuildCellGraph <- function(...) {
  id <- cli::cli_process_start("Computing PCA-based distance matrix ...")
  on.exit(cli::cli_process_done(id), add = TRUE)
  BuildCellGraph_inter(...)
}