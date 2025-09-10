
# ScGRADES: Single-cell Grading and Dataset Evaluation System

---

## Installation

### Create Conda environment

```bash
# Create a conda environment with required R packages
conda create -n scgrades_env -c conda-forge -c bioconda r-base=4.1 r-seurat r-dplyr r-tibble r-purrr r-optparse r-ggplot2
conda activate scgrades_env
```

### Install bigdist from GitHub (tested on v0.1.5)

```r
install.packages("remotes") # if not already installed
remotes::install_github("talegari/bigdist")
```

### Install scGRADES from CRAN

```r
#Install R package
install.packages("ScGRADES")
```

### Install scGRADES from GitHub

```r
# Install devtools if not already
install.packages("devtools")

# Install directly from GitHub
devtools::install_github("LibinLab/scGRADES")

#Install R package
install.packages("ScGRADES")
```

### Dependencies


```r
# R version:
R (4.1-4.2)

# R packages:
Seurat (≥ 4.0), bigdist(0.1.5), dplyr, purrr, tibble, optparse, ggplot2, cli, stats
```

---

## Example Data  

We provide a small Seurat object as an example dataset for testing and demonstration purposes.  This dataset has undergone **quality control (QC)**, **batch integration**, and has been processed with  `RunPCA()`. Therefore, it can be directly used for **ScGRADES**.  

**Download link:**  [example_scGRADRS.rds](https://drive.google.com/drive/u/1/folders/1JKjIb_qOeCEtV2kgGgAokmi4FgNACIdZ) 

---

## Usage 

```r
library(ScGRADES)

# Load the example dataset
seurat_obj <- readRDS("your_data.rds")
DefaultAssay(seurat_obj) <- "integrated"

# Step 1: Cluster the cells at multiple resolutions
seurat_obj <- FindClusterAcrossRes(seurat_obj, resolutions = seq(0.1, 1.5, 0.1))

# Step 2: Compute PCA-based Euclidean distance
dist_mat <- BuildCellGraph(seurat_obj, reduction = "pca")

# Step 3: Identify core cells
Seurat_obj <- IdentifyCoreCells(seurat_obj, dist = dist_mat, top_n = 20, resolutions = seq(0.1, 1.5, 0.1))
```

## Usage (Harmony integration)

```r
# Step 1: Cluster the cells at multiple resolutions
seurat_obj <- FindClusterAcrossRes(seurat_obj,resolutions = seq(0.1, 1.5, 0.1), reduction = "harmony")

# Step 2: Compute PCA-based Euclidean distance
dist_mat <- BuildCellGraph(seurat_obj, reduction = "harmony")

# Step 3: Identify core cells
Seurat_obj <- IdentifyCoreCells(seurat_obj, dist = dist_mat, top_n = 20, resolutions = seq(0.1, 1.5, 0.1))
```

---

## Output

- Adds a `CellPopulation` column to `seurat_obj@meta.data`, labeling cells as:
  - `core cell`: clusters are clear and consistent
  - `intermediate cell`: clusters are partially clear
  - `marginal cell`: clusters are unclear


---

## Citation
