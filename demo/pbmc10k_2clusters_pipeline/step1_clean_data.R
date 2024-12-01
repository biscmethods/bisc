#!/usr/bin/Rscript

library(here)  # To work with paths
# library(patchwork)
library(Seurat)
library(scregclust)
library(data.table)

# To be able to run library("EnsDb.Hsapiens.v79") you need to at least:
# BiocManager::install("GenomeInfoDb")
# BiocManager::install("SparseArray")
# BiocManager::install("EnsDb.Hsapiens.v79")
library("EnsDb.Hsapiens.v79")
# For sctransform to use speedy algos
# BiocManager::install('GenomicRanges', force=TRUE)
# BiocManager::install('glmGamPoi', force=TRUE)
library("glmGamPoi")

while (sink.number() > 0) {
  sink()
  sink(file = NULL)
}  # Because some of our scripts redirects output from scregclust to force it to be quite. This restores output.
gc()  # Force clean memory

# options(warn=2)  # To convert warning messages into error messages which display row of error. For debugging.

# Get absolute path where script is located, by using relative paths.
# Data from https://cellxgene.cziscience.com/collections/999f2a15-3d7e-440b-96ae-2c806799c08c
path_dataset <- file.path(local_data, "pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")
path_env_data <- file.path(local_data, "env_data_pbmc10k_sctransform.RData")
path_gex <- file.path(local_data, "GEX Graph-Based10k.csv")  # This you need to export from Loupe browser

if(redo_flag){
  if (file.exists(path_env_data)) {
    #Delete file if it exists
    file.remove(path_env_data)
  }
}


# Set seed for example
set.seed(214)

if (!file.exists(path_dataset)) {
  url <- "https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5"
  download.file(url, path_data, cacheOK = FALSE, mode = "wb")
}

if (file.exists(path_env_data)) {
  load(path_env_data)
}else {

  # Read data
  d <- Seurat::Read10X_h5(path_dataset, use.names = TRUE, unique.features = TRUE)[["Gene Expression"]]

  # Read cell cluster data from Loupe browser
  true_cluster_allocation <- read.csv(path_gex, header = TRUE)
  # Reorder true_cluster_allocation based on the column names of d
  true_cluster_allocation <- true_cluster_allocation[match(colnames(d), true_cluster_allocation$Barcode), ]
  true_cluster_allocation <- factor(as.numeric(gsub("Cluster ", "", true_cluster_allocation$GEX.Graph.based)))
  # Remove everything but cluster 1 and 2
  d <- d[,!is.na(true_cluster_allocation)]
  true_cluster_allocation <- true_cluster_allocation[!is.na(true_cluster_allocation)]



  d <- Seurat::CreateSeuratObject(counts = d, min.cells = 3, min.features = 200)

  # Filter
  d[["percent.mt"]] <- PercentageFeatureSet(d, pattern = "^MT.")
  d <- subset(d, subset = percent.mt < 30 & nFeature_RNA < 6000)

  # Read cell cluster data from Loupe browser
  true_cluster_allocation <- read.csv(path_gex, header = TRUE)
  # Reorder true_cluster_allocation based on the column names of d
  true_cluster_allocation <- true_cluster_allocation[match(colnames(d), true_cluster_allocation$Barcode), ]
  true_cluster_allocation <- factor(as.numeric(gsub("Cluster ", "", true_cluster_allocation$GEX.Graph.based)))



  # Pre-process
  options(future.globals.maxSize = 20000 * 1024^2)  # Standard was to only allow output of 500MiB. This is 3GiB
  future::plan("multicore", workers = 24)  # To set how many cores SCTransform will use?
  d <- Seurat::SCTransform(d, variable.features.n = 200)

  d <- d@assays$SCT$scale.data  # rows are genes, cols are cells
  # counts <- Seurat::LayerData(d, assay = "SCT", layer = "counts")
  # # Subset to same genes as after SCTransform normalization
  # counts <- counts[rownames(d), ]

  # Extract TFs and genesymbols
  out <- scregclust_format(d)

  genesymbols <- out$genesymbols
  is_regulator <- out$is_regulator
  rm(out)

  # pbmc_pearson_residuals <- cbind(data.frame(genesymbol = genesymbols), as.data.frame(expression))
  # pbmc_counts <- cbind(data.frame(genesymbol = rownames(counts)), as.data.frame(counts))
  # pbmc_genes_meta <- data.frame(genesymbol = genesymbols, is_regulator = is_regulator)




  # # Keep the correct type of cells
  # true_cluster_allocation <- d@meta.data$cell_type
  # keep_cells <- true_cluster_allocation == 'oligodendrocyte' | true_cluster_allocation == 'mature T cell'
  # keep_cells_split <- sample(c(1, 2), length(keep_cells), prob = c(0.15, 0.85), replace = T)
  # keep_ind <- keep_cells_split == 1
  # keep_cells <- keep_cells & keep_ind

  # true_cluster_allocation <- factor(true_cluster_allocation[keep_cells])
  # d <- d[, keep_cells]
  # rm(keep_cells)

  # Remove genes which are 'constant', aka no variance, just one number. Scregclust doesn't work on those.
  # non_constant_ind <- which(apply(d, 1, sd) != 0)
  # d <- d[non_constant_ind,]
  # rm(non_constant_ind)


  # Remove all genes with less non zero numbers than 200
  # d <- d[, apply(d, MARGIN = 2, function(x) sum(x != 0)) > 200]

  # Remove all genes/rows that don't correlate with other rows more than 0.1
  # cor_matrix <- abs(cor(t(d)))  # For correlation we want cols as genes
  # diag(cor_matrix) <- 0
  # threshold <- 0.1
  # keep_rows <- apply(cor_matrix, 1, max) > threshold
  # d <- d[keep_rows,]
  # rm(keep_rows, cor_matrix)

  # Ensembl gene IDs (like "ENSG00000269696") to standard gene names (like "A1BG"),
  # geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= rownames(d), keytype = "GENEID", columns = c("SYMBOL","GENEID"))
  # select_these <- geneIDs1[,2] # Since we probably lost genes in the translation
  # new_names <- geneIDs1[,1]
  # d <- d[select_these,]
  # rownames(d) <- new_names
  # rm(select_these)

  # Find out which genes/rows are regulators.
  # is_regulator <- give_me_regulators(normal_gene_names = rownames(d))

  # Put gene/rows in order. This puts the regulator genes at the end.
  d <- d[order(is_regulator),]
  is_regulator <- is_regulator[order(is_regulator)]

  # Put cells/columns in order.
  d <- d[, order(true_cluster_allocation)]
  true_cluster_allocation <- true_cluster_allocation[order(true_cluster_allocation)]

  n_regulator_genes <- sum(is_regulator)
  n_target_genes <- nrow(d) - n_regulator_genes

  unique_cell_types <- as.numeric(sort(unique(true_cluster_allocation)))
  n_cell_clusters <- length(unique_cell_types)
  # Now when we have put the cells in order we just need to count the cells
  # in each cell cluster. Then it's easy to make a vector with the true_cluster_allocation
  # further down in code
  n_cells_in_each_cluster <- numeric(n_cell_clusters)
  for(i_current_cell_cluster in seq(n_cell_clusters)){
    n_cells_in_each_cluster[i_current_cell_cluster] <- sum(true_cluster_allocation == unique_cell_types[i_current_cell_cluster])
  }
  n_total_cells <- sum(n_cells_in_each_cluster)

  # Print some stats
  print(object.size(d), units = 'MB', standard = 'SI')
  print(paste("Cells in clusters:"), quote = FALSE)
  print(paste("", unique_cell_types), quote = FALSE)
  print(paste("", n_cells_in_each_cluster), quote = FALSE)
  print(paste("Number of regulator genes are", sum(is_regulator)), quote = FALSE)
  print("Scregclust wants more cells than regulator genes x 2 for each cell cluster. Otherwise it doesn't work.", quote = FALSE)

  # Set variables ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # We assume cells are ordered in the order of cell clusters. So the first x columns are cell cluster 1, etc.
  true_cluster_allocation <- rep(1:n_cell_clusters, times = n_cells_in_each_cluster)
  rm(n_cells_in_each_cluster, unique_cell_types)

  ind_targetgenes <- which(c(rep(1, n_target_genes), rep(0, n_regulator_genes)) == 1)
  ind_reggenes <- which(c(rep(0, n_target_genes), rep(1, n_regulator_genes)) == 1)

  gc()  # Force clean memory

  save(d,
       is_regulator,
       n_cell_clusters,
       n_regulator_genes,
       n_target_genes,
       n_total_cells,
       ind_reggenes,
       ind_targetgenes,
       true_cluster_allocation,
       file = path_env_data)
}

