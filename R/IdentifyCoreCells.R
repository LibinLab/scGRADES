#' Identify Core Cells Based on Cluster Consistency
#'
#' @title Identify Core Cells
#' @description This function identifies "core" cells based on their local neighborhood consistency
#' across multiple clustering resolutions.
#' @name IdentifyCoreCells
#' @param scRNA A Seurat object.
#' @param dist A distance matrix (e.g., Euclidean distance between cells).
#' @param top_n Number of nearest neighbors to consider (default = 20).
#' @param resolutions A numeric vector of resolution values used in clustering (e.g., c(0.1, 0.2, ..., 1.5)).
#'
#' @return A Seurat object with updated meta.data including CellClass annotations.
#' @export
#' 
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate case_when select bind_rows
#' @importFrom tibble tibble
#' @importFrom rlang sym
#' @importFrom Seurat AddMetaData 
#' @importFrom bigdist bigdist_extract
#' @importFrom cli cli_process_start cli_process_done
#' 
#' @examples
#' \dontrun{
#' scRNA <- IdentifyCoreCells(scRNA, dist, top_n = 20, resolutions = seq(0.1, 1.5, 0.1))
#' }

utils::globalVariables(c("Cell1", "Cell2", "Dist", "CellClass", "pca"))

IdentifyCoreCells_inter <- function(scRNA, dist, top_n = 20, resolutions = seq(0.1, 1.5, 0.1)) {
#  library(dplyr)
#  library(purrr)
#  library(stringr)
  
  meta <- scRNA@meta.data
  pca <- Seurat::Embeddings(scRNA, reduction = "pca")
  cell_names <- rownames(pca)
  
  res_cols <- grep("snn_res", colnames(meta), value = TRUE)
  if (length(res_cols) == 0) stop("No clustering resolutions found (e.g., 'snn_res.') in metadata.")
  res_prefix <- stringr::str_extract(res_cols[1], "^.*snn_res\\.")
  
  # Get top-n neighbors for each cell
  top_list <- purrr::map(seq_along(cell_names), function(i) {
    #d_col <- dist[, i]
    d_col <- bigdist::bigdist_extract(dist, , i) 
    rownames(d_col) <- cell_names
    cell1 <- cell_names[i]
    df <- tibble::tibble(Cell2 = rownames(d_col), Dist = d_col[,1]) %>%
      dplyr::filter(Cell2 != cell1) %>%
      dplyr::slice_min(Dist, n = top_n)
    df$Cell1 <- cell1
    df
  })
  
  # Collect same-cluster ratio info across all resolutions
  all_ratios <- list()
  all_topn_annot <- list()
  
  for (res in resolutions) {
    #res_col <- paste0("integrated_snn_res.", res)
    res_col <- paste0(res_prefix, res)
    meta_res <- meta %>%
      dplyr::mutate(Cell1 = rownames(meta)) %>%
      dplyr::select(Cell1, Cluster = !!sym(res_col))
    # CellRatioFunc: how many top-n neighbors are in same cluster
    ratio_list <- purrr::map(top_list, function(df) {
      Cell1 <- df$Cell1[1]
      cl <- meta_res$Cluster[meta_res$Cell1 == Cell1]
      same_cluster_cells <- meta_res$Cell1[meta_res$Cluster == cl]
      n_match <- sum(df$Cell2 %in% same_cluster_cells)
      tibble::tibble(Cell1 = Cell1,
                     TopNSameClusterN = n_match,
                     Ratio = n_match / top_n)
    })
    annot_list <- purrr::map(top_list, function(df) {
      Cell1 <- df$Cell1[1]
      cl <- meta_res$Cluster[meta_res$Cell1 == Cell1]
      same_cluster_cells <- meta_res$Cell1[meta_res$Cluster == cl]
      df <- df %>% dplyr::mutate(SameClusterTF = Cell2 %in% same_cluster_cells)
      df$Resolution <- paste0("res_", res)
      df
    })
    res_ratio_df <- dplyr::bind_rows(ratio_list) %>%
      dplyr::mutate(Resolution = paste0("res_", res)) %>%
      dplyr::mutate(CellClass = case_when(
        Ratio == 1 ~ "Core",
        Ratio >= 0.5 ~ "Intermediate",
        TRUE ~ "Marginal"
      ))
    all_ratios[[as.character(res)]] <- res_ratio_df
    all_topn_annot[[as.character(res)]] <- dplyr::bind_rows(annot_list)
  }
  
  SameClusterRatio <- dplyr::bind_rows(all_ratios)
  SameClusterTopN <- dplyr::bind_rows(all_topn_annot)
  
  # Summarize CellClass across all resolutions
  sum_SameClusterRatio <- SameClusterRatio %>%
    dplyr::group_by(Cell1) %>%
    dplyr::summarise(
      SumCenter = sum(CellClass == "Core"),
      SumNormal = sum(CellClass == "Intermediate"),
      SumBorder = sum(CellClass == "Marginal")
    ) %>%
    dplyr::mutate(CellClass = dplyr::case_when(
      SumCenter == length(resolutions) ~ "core cell",
      SumBorder != 0 ~ "marginal cell",
      TRUE ~ "intermediate cell"
    ))
  # Add back to Seurat meta data
  CellPopulations <- sum_SameClusterRatio %>% 
    dplyr::select(Cell1, CellClass) %>% 
    dplyr::rename(CellPopulation = CellClass)
  
  CellPopulations <- as.data.frame(CellPopulations)
  rownames(CellPopulations) <- CellPopulations$Cell1
  CellPopulations <- CellPopulations %>% dplyr::select(-Cell1)
  scRNA <- Seurat::AddMetaData(scRNA, metadata = CellPopulations)
  return(scRNA)
}

IdentifyCoreCells <- function(...) {
  id <- cli::cli_process_start("Identify Core Cell population ...")
  on.exit(cli::cli_process_done(id), add = TRUE)
  IdentifyCoreCells_inter(...)
}