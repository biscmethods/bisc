#!/usr/bin/Rscript
rm(list = ls())

library(here)  # To work with paths
library(patchwork)
sink()

# options(warn=2)  # To convert warning messages into error messages which display row of error. For debugging.

# Get absolute path where script is located, by using relative paths.
demo_path <- here::here("demo")
R_path <- here::here("R")
output_path <- demo_path

source(file.path(R_path, "generate_dummy_data_for_cell_clustering.R"))
source(file.path(R_path, "biclust_scregclust.R"))
source(file.path(R_path, "randomise_cluster_labels.R"))

#############################################
############ data for dev ###################
#############################################

# Set seed for example
set.seed(214)

# Load data
# Folders and filenames
path_data <- here::here('data')
path_Neftel2019 <- file.path(path_data, "Neftel2019")
path_Group2 <- file.path(path_Neftel2019, "Group2")
path_neftel_seurat_group2 <- file.path(path_data, "neftel_seurat_group2.rds")

# Load Neftel 2019 data prepped in download_and_prep_data.R "group 2"
neftel_smartseq2 <- readRDS(path_neftel_seurat_group2)
neftel_smartseq2 <- neftel_smartseq2@assays$SCT@scale.data

# Remove all genes/rows that don't correlate with other rows more than 0.1
cor_matrix <- abs(cor(t(d)))
diag(cor_matrix) <- 0
threshold <- 0.1
keep_rows <- apply(cor_matrix, 1, max) > threshold
neftel_smartseq2 <- neftel_smartseq2[keep_rows,]

# Remove all cells that aren't malignant
cells <- read.table(file = file.path(path_Group2, 'cells.csv'), sep = ' ', header = TRUE, stringsAsFactors = FALSE)
cells_malignant <- cells[cells[, 'malignant'] == 'yes',]
neftel_smartseq2_malignant <- neftel_smartseq2[, colnames(neftel_smartseq2) %in% cells_malignant[, 'cell_name']]

# Figure out which genes are regulator using the names
out <- scregclust::scregclust_format(neftel_smartseq2_malignant)
is_regulator <- out[['is_regulator']]

print(paste("Out of", ncol(neftel_smartseq2), "cells we selected", ncol(neftel_smartseq2_malignant), "that are malignant."), quote = FALSE)
print(paste("Out of", length(is_regulator), "genes", sum(is_regulator == 1), "are regulator genes."), quote = FALSE)
print(paste("The dimension of the matrix is", nrow(neftel_smartseq2_malignant), "rows and", ncol(neftel_smartseq2_malignant), "cols."), quote = FALSE)

# Put target genes first, then regulator genes
neftel_smartseq2_malignant_ordered <- neftel_smartseq2_malignant[order(is_regulator),]
is_regulator <- is_regulator[order(is_regulator)]

# True cell cluster allocation
# We need to order the columns by the cell clusters
# We will use MES, AC, OPC and NPC as true cell clusters
df <- cells_malignant[c("cell_name", "MESlike2", "MESlike1", "AClike", "OPClike", "NPClike1", "NPClike2")]
df <- df[complete.cases(df),]
colors <- apply(df[, 2:7], MARGIN = c(1), which.max)
df['cell_cluster'] <- c("MES", "MES", "AC", "OPC", "NPC", "NPC")[colors]

# Do a plot if you want
# pca_res <- prcomp(df[, 2:7], scale. = TRUE)
# ggfortify::autoplot(pca_res, data = df, colour = 'cell_cluster')

print(paste("Number of duplicated cell names in df is", sum(duplicated(df$cell_name))), quote = FALSE)
df <- df[!duplicated(df$cell_name),]

# There are more cells/columns in neftel gene data matrix then there is in
# the cells.csv file. This will fix remove the mis-match.
a <- df$cell_name
b <- colnames(neftel_smartseq2_malignant_ordered)
neftel_smartseq2_malignant_exists_in_cells <- neftel_smartseq2_malignant_ordered[, b %in% a]

# This is the magic line that actual puts the cells/columns
# in the order of cell clusters
df <- df[order(df$cell_cluster),]

# Then we need to apply that order to the gene matrix
neftel_smartseq2_malignant_sorted <- neftel_smartseq2_malignant_exists_in_cells[, df$cell_name]

# Now when we have put the cells in order we just need to count the cells
# in each cell cluster. Then it's easy to make a vector with the true_cluster_allocation
# further down in code
n_cells_cell_cluster_1 <- sum(df$cell_cluster == 'AC')
n_cells_cell_cluster_2 <- sum(df$cell_cluster == 'MES')
n_cells_cell_cluster_3 <- sum(df$cell_cluster == 'NPC')
n_cells_cell_cluster_4 <- sum(df$cell_cluster == 'OPC')

print(paste("Cells in clusters:"), quote = FALSE)
print(paste("1:", n_cells_cell_cluster_1), quote = FALSE)
print(paste("2:", n_cells_cell_cluster_2), quote = FALSE)
print(paste("3:", n_cells_cell_cluster_3), quote = FALSE)
print(paste("4:", n_cells_cell_cluster_4), quote = FALSE)
print(paste("Number of regulator genes are", sum(is_regulator)), quote = FALSE)
print("Scregclust wants more cells than regulator genes x 2 for each cell cluster. Otherwise it doesn't work.", quote = FALSE)

# Double check the sums are the same
total_n_cells_in_cellstxt <- n_cells_cell_cluster_1 +
  n_cells_cell_cluster_2 +
  n_cells_cell_cluster_3 +
  n_cells_cell_cluster_4
total_n_cells_in_genedata <- ncol(neftel_smartseq2_malignant_sorted)
print(paste(str(total_n_cells_in_cellstxt), " = ", str(total_n_cells_in_genedata)), quote = FALSE)


#--------------------------


# Set variables ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

n_target_genes <- length(is_regulator) - sum(is_regulator == 1)
n_regulator_genes <- sum(is_regulator == 1)
n_total_cells <- ncol(neftel_smartseq2_malignant_sorted)

# We assume cells are ordered in the order of cell clusters. So the first x columns are cell cluster 1, etc.
n_cells <- c(n_cells_cell_cluster_1, n_cells_cell_cluster_2, n_cells_cell_cluster_3, n_cells_cell_cluster_4)
true_cluster_allocation <- rep(1:n_cell_clusters, times = n_cells)  # TODO: do this

# Because "dat <- cbind(Z_t, Z_r)" in generate_dummy_data_for_cell_clustering
ind_targetgenes <- which(c(rep(1, n_target_genes), rep(0, n_regulator_genes)) == 1)
ind_reggenes <- which(c(rep(0, n_target_genes), rep(1, n_regulator_genes)) == 1)

# Set up some variables
n_cell_clusters <- length(unique(disturbed_initial_cell_clust))
n_target_genes <- length(ind_targetgenes)
n_regulator_genes <- length(ind_reggenes)

cell_id <- 1:n_total_cells


min_number_of_clusters <- 10
max_number_of_clusters <- 10
target_gene_cluster_vector <- seq(min_number_of_clusters, max_number_of_clusters)
gc()


# Run screg with a bunch of different cluster setings
results <- vector(mode = "list", length = length(target_gene_cluster_vector))

current_cell_cluster <- neftel_smartseq2_malignant_sorted[, df$cell_cluster == 'AC']
for (i_n_target_genes_clusters in seq(length(target_gene_cluster_vector))) {
  n_target_genes_clusters <- target_gene_cluster_vector[i_n_target_genes_clusters]
  print(paste("Cell cluster", i_cell_cluster, "Number of target gene clusters", n_target_genes_clusters), quote = FALSE)

  results[[i_n_target_genes_clusters]] <- scregclust::scregclust(
    expression = current_cell_cluster,  # p rows of genes, n columns of cells
    split_indices = NULL,
    genesymbols = 1:(n_target_genes + n_regulator_genes),
    is_regulator = is_regulator,
    n_cl = n_target_genes_clusters,
    penalization = 0.0001,
    noise_threshold = 0.000001,
    verbose = TRUE,
    n_cycles = 40,
    compute_silhouette = FALSE,
    center = FALSE
  )
}
