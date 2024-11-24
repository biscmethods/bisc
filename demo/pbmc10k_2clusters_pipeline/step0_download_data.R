#!/usr/bin/Rscript
rm(list = ls())

library(here)  # To work with paths
library(Seurat)
library(scregclust)
library(data.table)
# BiocManager::install('GenomicRanges', force=TRUE)
library("glmGamPoi")

# Get absolute path where script is located, by using relative paths.
demo_path <- here::here("demo")
R_path <- here::here("R")
output_path <- demo_path
path_data <- here::here('data')

path_dataset <- file.path(path_data, "pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")
path_env_data <- file.path(path_data, "env_data_2clusters_sctransform.RData")
path_general_env_data <- file.path(path_data, "env_data_general_2clusters_sctransform.RData")
path_cluster1 <- file.path(path_data, "env_data_cluster1.RData")
path_cluster2 <- file.path(path_data, "env_data_cluster2.RData")


if (!file.exists(path_dataset)) {
  url <- "https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5"
  download.file(url, path_data, cacheOK = FALSE, mode = "wb")
}

pbmc_data <- Read10X_h5(
  path_dataset, use.names = TRUE, unique.features = TRUE
)[["Gene Expression"]]

pbmc <- CreateSeuratObject(
  counts = pbmc_data, min.cells = 3, min.features = 200
)

# Filter
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT.")
pbmc <- subset(pbmc, subset = percent.mt < 30 & nFeature_RNA < 6000)
# Pre-process
options(future.globals.maxSize = 20000 * 1024^2)  # Standard was to only allow output of 500MiB. This is 3GiB
pbmc <- SCTransform(pbmc, variable.features.n = 2000)

expression <- LayerData(pbmc, assay = "SCT", layer = "scale.data")
counts <- LayerData(pbmc, assay = "SCT", layer = "counts")
# Subset to same genes as after SCTransform normalization
counts <- counts[rownames(expression), ]

# Extract TFs and genesymbols
out <- scregclust_format(expression)

genesymbols <- out$genesymbols
is_regulator <- out$is_regulator

fwrite(
  cbind(
    data.frame(genesymbol = rownames(expression)),
    as.data.frame(expression)
  ),
  file = "data/pbmc-pearson-resid.csv",
  quote = FALSE
)
fwrite(
  cbind(
    data.frame(genesymbol = rownames(counts)),
    as.data.frame(counts)
  ),
  file = "data/pbmc-counts.csv",
  quote = FALSE
)
fwrite(
  data.frame(
    genesymbol = genesymbols,
    is_regulator = is_regulator
  ),
  file = "data/pbmc-genes-meta.csv",
  quote = FALSE
)
