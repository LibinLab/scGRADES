#' Identify Core Cells Across Clustering Resolutions
#'
#' Identify core cells from K-nearest-neighbor consistency across multiple
#' clustering resolutions. At a given resolution, a cell is considered core
#' when all of its K nearest neighbors belong to the same cluster. The final
#' classification requires the fraction of resolutions at which the cell is
#' core to be greater than or equal to a user-specified CS cutoff.
#'
#' @param scRNA A Seurat object containing clustering results in `meta.data`.
#' @param dist A `bigdist` distance object whose cell order is identical to the
#'   row order of `scRNA@meta.data`.
#' @param top_n A positive integer specifying the number of nearest neighbors,
#'   K. The default is `20`.
#' @param resolutions A non-empty numeric vector containing the clustering
#'   resolutions to evaluate. The default is `seq(0.1, 1.5, 0.1)`.
#' @param CS One or more numeric cutoffs between `0` and `1`. A cell is labeled
#'   as a core cell for a cutoff when the fraction of resolutions at which it
#'   is core is greater than or equal to that cutoff. The default `CS = 1`
#'   reproduces the strict manuscript definition, which requires the cell to
#'   be core at every resolution.
#'
#' @return The input Seurat object with one `CellPopulation_CS_<cutoff>` column
#'   added to `meta.data` for each requested CS cutoff. Each column contains
#'   either `"core cell"` or `"marginal cell"`. For example, `CS = 0.95`
#'   produces `CellPopulation_CS_0.95`.
#'
#' @details
#' For cell i at resolution r, the within-cluster neighbor consistency is
#'
#' `S(r, i) = (number of same-cluster neighbors) / K`.
#'
#' The cell is core at resolution r only when `S(r, i) = 1`. For a user-supplied
#' cutoff, the final classification is core when
#'
#' `(number of core resolutions) / (number of evaluated resolutions) >= CS`.
#'
#' Exactly `top_n` neighbors are retained. If distances are tied at the K-th
#' position, ties are resolved using their existing cell order.
#'
#' @name IdentifyCoreCells
#'
#' @importFrom bigdist bigdist_extract
#' @importFrom cli cli_process_done cli_process_start
#' @importFrom dplyr filter slice_min
#' @importFrom Seurat AddMetaData
#' @importFrom stringr str_extract
#' @importFrom tibble tibble
#'
#' @examples
#' \dontrun{
#' # Strict manuscript definition
#' seurat_obj <- IdentifyCoreCells(
#'   scRNA = seurat_obj,
#'   dist = dist_mat,
#'   top_n = 20,
#'   resolutions = seq(0.1, 1.5, 0.1),
#'   CS = 1
#' )
#'
#' # Robustness analysis using multiple CS cutoffs
#' seurat_obj <- IdentifyCoreCells(
#'   scRNA = seurat_obj,
#'   dist = dist_mat,
#'   top_n = 20,
#'   resolutions = seq(0.1, 1.5, 0.1),
#'   CS = c(1, 0.95, 0.9, 0.8)
#' )
#' }
NULL

utils::globalVariables(c("Cell2", "Dist"))

# Internal implementation used by the exported wrapper.
IdentifyCoreCells_inter <- function(scRNA,
                                    dist,
                                    top_n = 20,
                                    resolutions = seq(0.1, 1.5, 0.1),
                                    CS = 1) {
  if (!inherits(scRNA, "Seurat")) {
    stop("scRNA must be a Seurat object.", call. = FALSE)
  }

  if (!is.numeric(CS) ||
      length(CS) < 1L ||
      anyNA(CS) ||
      any(!is.finite(CS)) ||
      any(CS < 0 | CS > 1)) {
    stop(
      "CS must contain numeric values between 0 and 1.",
      call. = FALSE
    )
  }
  CS <- unique(as.numeric(CS))

  if (!is.numeric(top_n) ||
      length(top_n) != 1L ||
      is.na(top_n) ||
      !is.finite(top_n) ||
      top_n < 1 ||
      top_n != as.integer(top_n)) {
    stop("top_n must be one positive integer.", call. = FALSE)
  }
  top_n <- as.integer(top_n)

  if (!is.numeric(resolutions) ||
      length(resolutions) < 1L ||
      anyNA(resolutions) ||
      any(!is.finite(resolutions))) {
    stop(
      "resolutions must be a non-empty numeric vector without missing values.",
      call. = FALSE
    )
  }
  resolutions <- unique(as.numeric(resolutions))

  meta <- scRNA@meta.data
  cell_names <- rownames(meta)

  if (length(cell_names) < 2L || anyDuplicated(cell_names)) {
    stop(
      "scRNA must contain at least two cells with unique cell names.",
      call. = FALSE
    )
  }

  if (top_n >= length(cell_names)) {
    stop(
      "top_n must be smaller than the number of cells.",
      call. = FALSE
    )
  }

  res_cols <- grep("snn_res", colnames(meta), value = TRUE)
  if (length(res_cols) == 0L) {
    stop(
      "No clustering resolution columns were found in scRNA@meta.data.",
      call. = FALSE
    )
  }

  res_prefix <- stringr::str_extract(res_cols[1], "^.*snn_res\\.")
  if (is.na(res_prefix)) {
    stop(
      "Unable to determine the clustering-column prefix from metadata.",
      call. = FALSE
    )
  }

  expected_res_cols <- paste0(res_prefix, resolutions)
  missing_res_cols <- setdiff(expected_res_cols, colnames(meta))
  if (length(missing_res_cols) > 0L) {
    stop(
      paste(
        "Missing clustering columns:",
        paste(missing_res_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (anyNA(meta[, expected_res_cols, drop = FALSE])) {
    stop(
      "Clustering resolution columns must not contain missing labels.",
      call. = FALSE
    )
  }

  top_list <- lapply(seq_along(cell_names), function(i) {
    d_col <- bigdist::bigdist_extract(dist, , i)

    if (NROW(d_col) != length(cell_names)) {
      stop(
        "dist and scRNA contain different numbers of cells.",
        call. = FALSE
      )
    }

    rownames(d_col) <- cell_names
    cell1 <- cell_names[i]

    df <- tibble::tibble(
      Cell2 = rownames(d_col),
      Dist = d_col[, 1]
    )
    df <- dplyr::filter(df, Cell2 != cell1)
    df <- dplyr::slice_min(
      df,
      Dist,
      n = top_n,
      with_ties = FALSE
    )
    df$Cell1 <- cell1
    df
  })

  all_ratios <- vector("list", length(resolutions))

  for (res_index in seq_along(resolutions)) {
    res <- resolutions[res_index]
    res_col <- paste0(res_prefix, res)

    meta_res <- data.frame(
      Cell1 = cell_names,
      Cluster = meta[[res_col]],
      stringsAsFactors = FALSE
    )

    ratio_list <- lapply(top_list, function(df) {
      cell1 <- df$Cell1[1]
      cluster_i <- meta_res$Cluster[meta_res$Cell1 == cell1]
      same_cluster_cells <- meta_res$Cell1[
        meta_res$Cluster == cluster_i
      ]
      n_match <- sum(df$Cell2 %in% same_cluster_cells)

      data.frame(
        Cell1 = cell1,
        TopNSameClusterN = n_match,
        Ratio = n_match / top_n,
        Resolution = paste0("res_", res),
        stringsAsFactors = FALSE
      )
    })

    all_ratios[[res_index]] <- do.call(rbind, ratio_list)
  }

  same_cluster_ratio <- do.call(rbind, all_ratios)

  core_resolution_n <- tapply(
    same_cluster_ratio$TopNSameClusterN == top_n,
    same_cluster_ratio$Cell1,
    sum
  )
  resolution_n <- table(same_cluster_ratio$Cell1)

  core_resolution_n <- as.integer(core_resolution_n[cell_names])
  resolution_n <- as.integer(resolution_n[cell_names])

  if (anyNA(core_resolution_n) ||
      anyNA(resolution_n) ||
      any(resolution_n != length(resolutions))) {
    stop(
      "Some cells do not have results at every resolution.",
      call. = FALSE
    )
  }

  metadata <- data.frame(
    .row_id = seq_along(cell_names),
    row.names = cell_names,
    stringsAsFactors = FALSE
  )
  metadata$.row_id <- NULL

  format_cs_label <- function(x) {
    label <- formatC(x, format = "f", digits = 6)
    sub("\\.?0+$", "", label)
  }

  cs_labels <- vapply(CS, format_cs_label, character(1))
  if (anyDuplicated(cs_labels)) {
    stop(
      "CS values must remain distinct when rounded to six decimal places.",
      call. = FALSE
    )
  }

  for (cs_index in seq_along(CS)) {
    cs_cutoff <- CS[cs_index]
    required_resolutions <- ceiling(
      cs_cutoff * length(resolutions) - sqrt(.Machine$double.eps)
    )
    cs_result <- as.integer(
      core_resolution_n >= required_resolutions
    )
    population <- ifelse(
      cs_result == 1L,
      "core cell",
      "marginal cell"
    )

    output_col <- paste0(
      "CellPopulation_CS_",
      cs_labels[cs_index]
    )
    metadata[[output_col]] <- population
  }

  scRNA <- Seurat::AddMetaData(scRNA, metadata = metadata)
  return(scRNA)
}

#' @rdname IdentifyCoreCells
#' @export
IdentifyCoreCells <- function(scRNA,
                              dist,
                              top_n = 20,
                              resolutions = seq(0.1, 1.5, 0.1),
                              CS = 1) {
  id <- cli::cli_process_start("Identifying core-cell population ...")
  on.exit(cli::cli_process_done(id), add = TRUE)

  IdentifyCoreCells_inter(
    scRNA = scRNA,
    dist = dist,
    top_n = top_n,
    resolutions = resolutions,
    CS = CS
  )
}
