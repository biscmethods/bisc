#!/usr/bin/Rscript
rm(list = ls())

library(here)  # To work with paths
# library(patchwork)
library(scregclust)
# To be able to run library("EnsDb.Hsapiens.v79") you need to at least:
# BiocManager::install("GenomeInfoDb")
# BiocManager::install("SparseArray")
# BiocManager::install("EnsDb.Hsapiens.v79")
library("EnsDb.Hsapiens.v79")
sink()  # Because some of our scripts redirects output from scregclust to force it to be quite. This restores output.
gc()  # Force clean memory

# options(warn=2)  # To convert warning messages into error messages which display row of error. For debugging.

# Get absolute path where script is located, by using relative paths.
demo_path <- here::here("demo")
R_path <- here::here("R")
output_path <- demo_path

path_data <- here::here('data')

# Data from https://cellxgene.cziscience.com/collections/999f2a15-3d7e-440b-96ae-2c806799c08c
path_coreGBmap <- file.path(path_data, "medium.rds")
path_env_data <- file.path(path_data, "env_data_2clusters_sctransform.RData")
path_general_env_data <- file.path(path_data, "env_data_general_2clusters_sctransform.RData")
path_cluster1 <- file.path(path_data, "env_data_cluster1.RData")
path_cluster2 <- file.path(path_data, "env_data_cluster2.RData")


# Set seed for example
set.seed(214)

# Ensemble gene IDs (like "ENSG00000269696") to standard gene names (like "A1BG"),
ensembl_to_normal <- function(data_matrix, ensembl_gene_names) {
  geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys = ensembl_gene_names, keytype = "GENEID", columns = c("SYMBOL", "GENEID"))
  select_these <- geneIDs1[, 2] # Since we probably lost genes in the translation
  normal_gene_names <- geneIDs1[, 1]

  if (nrow(data_matrix) == ncol(data_matrix)) {
    stop("Can't determine which dimension is genes/cells.")
  }

  if (nrow(data_matrix) == length(ensembl_gene_names)) {
    data_matrix <- data_matrix[select_these,]
  }else if (ncol(data_matrix) == length(ensembl_gene_names)) {
    data_matrix <- data_matrix[, select_these]
  }else(
    stop("The input ensembl vector doesn't match any dimension of data_matrix")
  )

  return(list("data_matrix" = data_matrix, "normal_gene_names" = normal_gene_names))
}

give_me_regulators <- function(normal_gene_names) {
  fake_matrix <- matrix(0, nrow = length(normal_gene_names), ncol = 1)  # This is fast
  rownames(fake_matrix) <- normal_gene_names  # as.vector(d@assays$SCT@var.features)
  out <- scregclust::scregclust_format(fake_matrix)  # Needs to be a matrix to use this
  is_regulator <- out[['is_regulator']]
  rm(fake_matrix, out)
  return(is_regulator)
}

if (file.exists(path_env_data) && file.exists(path_general_env_data)) {
  load(path_env_data)
}else {
  # Read data
  d <- readRDS(path_coreGBmap)

  # Keep the correct type of cells
  cell_types <- d@meta.data$cell_type
  keep_cells <- cell_types == 'oligodendrocyte' | cell_types == 'mature T cell'
  keep_cells_split <- sample(c(1, 2), length(keep_cells), prob = c(0.15, 0.85), replace = T)
  keep_ind <- keep_cells_split == 1
  keep_cells <- keep_cells & keep_ind

  cell_types <- factor(cell_types[keep_cells])
  d <- d[, keep_cells]
  rm(keep_cells)

  # Remove genes which are 'constant', aka no variance, just one number. Scregclust doesn't work on those.
  non_constant_ind <- which(apply(d@assays$RNA$data, 1, sd) != 0)
  d <- d[non_constant_ind,]
  rm(non_constant_ind)

  # Reduce the number of genes with sctransform. It converts gene names to ensembl gene names.
  options(future.globals.maxSize = 20000 * 1024^2)  # Standard was to only allow output of 500MiB or something. This is 3GiB
  future::plan("multicore", workers = 24)  # To set how many cores SCTransform will use?
  d <- Seurat::SCTransform(d, variable.features.n = 1000)

  d <- d@assays$SCT$scale.data  # rows are genes, cols are cells

  # Remove all genes with less non zero numbers than 200
  d <- d[, apply(d, MARGIN = 2, function(x) sum(x != 0)) > 200]

  # Remove all genes/rows that don't correlate with other rows more than 0.1
  # cor_matrix <- abs(cor(t(d)))  # For correlation we want cols as genes
  # diag(cor_matrix) <- 0
  # threshold <- 0.1
  # keep_rows <- apply(cor_matrix, 1, max) > threshold
  # d <- d[keep_rows,]
  # rm(keep_rows, cor_matrix)

  # Ensembl gene IDs (like "ENSG00000269696") to standard gene names (like "A1BG"),
  geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= rownames(d), keytype = "GENEID", columns = c("SYMBOL","GENEID"))
  select_these <- geneIDs1[,2] # Since we probably lost genes in the translation
  new_names <- geneIDs1[,1]
  d <- d[select_these,]
  rownames(d) <- new_names
  rm(select_these)

  # Find out which genes/rows are regulators.
  is_regulator <- give_me_regulators(normal_gene_names = rownames(d))

  # Put gene/rows in order. This puts the regulator genes at the end.
  d <- d[order(is_regulator),]
  is_regulator <- is_regulator[order(is_regulator)]

  # Put cells/columns in order.
  d <- d[, order(cell_types)]
  cell_types <- cell_types[order(cell_types)]

  n_regulator_genes <- sum(is_regulator)
  n_target_genes <- nrow(d) - n_regulator_genes


  # Now when we have put the cells in order we just need to count the cells
  # in each cell cluster. Then it's easy to make a vector with the true_cluster_allocation
  # further down in code
  n_cells_cell_cluster_1 <- sum(cell_types == 'oligodendrocyte')
  n_cells_cell_cluster_2 <- sum(cell_types == 'mature T cell')

  # Print some stats
  print(object.size(d), units = 'MB', standard = 'SI')
  print(paste("Cells in clusters:"), quote = FALSE)
  print(paste("oligodendrocyte:", n_cells_cell_cluster_1), quote = FALSE)
  print(paste("mature T cell:", n_cells_cell_cluster_2), quote = FALSE)
  print(paste("Number of regulator genes are", sum(is_regulator)), quote = FALSE)
  print("Scregclust wants more cells than regulator genes x 2 for each cell cluster. Otherwise it doesn't work.", quote = FALSE)
  rm(cell_types)

  # Set variables ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  n_cell_clusters <- 2
  # We assume cells are ordered in the order of cell clusters. So the first x columns are cell cluster 1, etc.
  n_cells_in_each_cluster <- c(n_cells_cell_cluster_1, n_cells_cell_cluster_2)
  true_cluster_allocation <- rep(1:n_cell_clusters, times = n_cells_in_each_cluster)
  rm(n_cells_in_each_cluster, n_cells_cell_cluster_1, n_cells_cell_cluster_2)

  ind_targetgenes <- which(c(rep(1, n_target_genes), rep(0, n_regulator_genes)) == 1)
  ind_reggenes <- which(c(rep(0, n_target_genes), rep(1, n_regulator_genes)) == 1)

  gc()  # Force clean memory


  save(d,
       is_regulator,
       n_cell_clusters,
       n_regulator_genes,
       n_target_genes,
       ind_reggenes,
       ind_targetgenes,
       true_cluster_allocation,
       file = path_env_data)

   save(is_regulator,
       n_cell_clusters,
       n_regulator_genes,
       n_target_genes,
       ind_reggenes,
       ind_targetgenes,
       true_cluster_allocation,
       file = path_general_env_data)
}

# Create seperate datasets for each cell state
if (!file.exists(path_cluster1) ||
  !file.exists(path_cluster2)) {
  print(object.size(d), units = 'MB', standard = 'SI')

  load(path_env_data)

  # Save one variable per cell cluster to save on memory use
  cell_cluster_1 <- d[, true_cluster_allocation == 1]
  print(dim(cell_cluster_1))
  print(object.size(cell_cluster_1), units = 'MB', standard = 'SI')
  save(cell_cluster_1, file = path_cluster1)
  rm(cell_cluster_1)
  gc()  # Force clean memory

  cell_cluster_2 <- d[, true_cluster_allocation == 2]
  print(dim(cell_cluster_2))
  print(object.size(cell_cluster_2), units = 'MB', standard = 'SI')
  save(cell_cluster_2, file = path_cluster2)
  rm(cell_cluster_2)
  gc()  # Force clean memory
}
