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
path_neftel_mn_group2 <- file.path(path_data, "neftel_mn_group2.rds")

neftel_smartseq2 <- readRDS(path_neftel_mn_group2)
cells <- read.table(file = file.path(path_Group2, 'cells.csv'), sep = ' ', header = TRUE, stringsAsFactors = FALSE)
cells_malignant <- cells[cells[, 'malignant'] == 'yes',]
neftel_smartseq2_malignant <- neftel_smartseq2[, colnames(neftel_smartseq2) %in% cells_malignant[, 'cell_name']]
out <- scregclust::scregclust_format(neftel_smartseq2_malignant)
is_regulator <- out[['is_regulator']]

print(paste("Out of", ncol(neftel_smartseq2), "cells we selected", ncol(neftel_smartseq2_malignant), "that are malignant."), quote = FALSE)
print(paste("Out of", length(is_regulator), "genes", sum(is_regulator == 1), "are regulator genes."), quote = FALSE)
print(paste("The dimension of the matrix is", nrow(neftel_smartseq2_malignant), "rows and", ncol(neftel_smartseq2_malignant), "cols."), quote = FALSE)

# Put some stuff in order
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


# Double check the sums are the same
total_n_cells_in_cellstxt <- n_cells_cell_cluster_1 +
  n_cells_cell_cluster_2 +
  n_cells_cell_cluster_3 +
  n_cells_cell_cluster_4
total_n_cells_in_genedata <- ncol(neftel_smartseq2_malignant_sorted)
print(paste(str(total_n_cells_in_cellstxt), " = ", str(total_n_cells_in_genedata)), quote = FALSE)


# Set variables ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

n_cell_clusters <- 4
n_target_gene_clusters <- c(3, 3, 3, 3)  # TODO:Number of target gene clusters in each cell cluster
n_target_genes <- length(is_regulator) - sum(is_regulator == 1)
n_regulator_genes <- sum(is_regulator == 1)
n_total_cells <- ncol(neftel_smartseq2_malignant_sorted)
# We assume cells are ordered in the order of cell clusters. So the first x columns are cell cluster 1, etc.
n_cells <- c(n_cells_cell_cluster_1, n_cells_cell_cluster_2, n_cells_cell_cluster_3, n_cells_cell_cluster_4)
true_cluster_allocation <- rep(1:n_cell_clusters, times = n_cells)  # TODO: do this


# Because "dat <- cbind(Z_t, Z_r)" in generate_dummy_data_for_cell_clustering
ind_targetgenes <- which(c(rep(1, n_target_genes), rep(0, n_regulator_genes)) == 1)
ind_reggenes <- which(c(rep(0, n_target_genes), rep(1, n_regulator_genes)) == 1)

disturbed_initial_cell_clust <- factor(randomise_cluster_labels(cluster_labels = true_cluster_allocation,
                                                                fraction_randomised = 0.25))

biclust_input_data <- t(neftel_smartseq2_malignant_sorted)
colnames(biclust_input_data) <- c(paste0("t", 1:n_target_genes), paste0("r", 1:n_regulator_genes))
biclust_input_data <- tibble::as_tibble(biclust_input_data)

# # These needs to be strings for discrete labels in pca plot
# data_for_plotting <- tibble::as_tibble(true_cell_cluster_allocation = generated_data$true_cell_clust,
#                                        biclust_input_data)
# pca_res <- prcomp(biclust_input_data, scale. = TRUE)
# p <- ggplot2::autoplot(pca_res, data = data_for_plotting, colour = 'true_cell_cluster_allocation')

# Set up some variables
n_cell_clusters <- length(unique(disturbed_initial_cell_clust))
n_target_genes <- length(ind_targetgenes)
n_regulator_genes <- length(ind_reggenes)


#############################################
############ end data for dev ###############
#############################################

###########################################
####initialise variables for dev ##########
##########################################

max_iter <- 100
initial_clustering <- disturbed_initial_cell_clust
n_target_genes <- n_target_genes
n_regulator_genes <- n_regulator_genes
n_total_cells <- n_total_cells
n_cell_clusters <- n_cell_clusters
ind_targetgenes <- ind_targetgenes
ind_reggenes <- ind_reggenes
# output_path <- modded_output_path

i_cell_cluster <- 1

use_weights <- TRUE
use_complex_cluster_allocation <- FALSE

demo_path <- here::here("demo")
output_path <- demo_path
# initial_clustering[1:14900] <- rep(1, 14900)
# initial_clustering[14901:15000] <- rep(2, 100)
# print(length(initial_clustering))
# stop("hej")


###########################################
###END initialise variables for dev #######
###########################################
# Split data into train/test
cell_data_split <- sample(c(1, 2), nrow(biclust_input_data), prob = c(0.5, 0.5), replace = T)
train_indices <- which(cell_data_split == 1)
test_indices <- which(cell_data_split == 2)

biclust_input_data_train <- biclust_input_data[train_indices,]
biclust_input_data_test <- biclust_input_data[test_indices,]

# TODO: Put this inside generate_data_lm or something
# Setup variables that will be used throughout the script
# We assume target genes start with t then a number. And likewise for regulator genes.
n_total_cells_train <- length(train_indices)
cell_id <- 1:n_total_cells
cell_id_train <- cell_id[train_indices]

initial_clustering_train <- initial_clustering[train_indices]
if (length(unique(initial_clustering)) != length(unique(initial_clustering_train))) {
  stop("Training set doesn't have all the clusters")
}

# Set up some variables
n_cell_clusters_train <- length(unique(initial_clustering_train))

# true_cell_cluster_allication <- factor(generated_data$true_cell_clust)
# true_cell_cluster_allication_train <- true_cell_cluster_allication[train_indices]
# true_cell_cluster_allication_train <- true_cell_cluster_allication[train_indices]

penalization_lambdas <- c(0.00001, 0.001, 0.1, 0.3, 0.5, 0.7)
BICLUST_RESULTS <- vector(mode = "list", length = length(penalization_lambdas))


for (i_penalization_lambda in seq_along(penalization_lambdas)) {
  print("", quote = FALSE)
  print(paste("Running biclust for penalization_lambda", penalization_lambdas[i_penalization_lambda]), quote = FALSE)
  BICLUST_RESULTS[[i_penalization_lambda]] <- biclust(dat = biclust_input_data,
                                                      cell_id = cell_id,
                                                      true_cell_cluster_allocation = factor(true_cluster_allocation),
                                                      max_iter = max_iter,
                                                      n_target_gene_clusters,
                                                      initial_clustering,
                                                      n_target_genes,
                                                      n_regulator_genes,
                                                      n_total_cells,
                                                      n_cell_clusters,
                                                      ind_targetgenes,
                                                      ind_reggenes,
                                                      output_path,
                                                      penalization_lambda = penalization_lambdas[i_penalization_lambda])

}

print("", quote = FALSE)
print("", quote = FALSE)
for (i_penalization_lambda in seq_along(penalization_lambdas)) {
  if (is.na(BICLUST_RESULTS[i_penalization_lambda])) {
    print(paste("penalization_lambda", penalization_lambdas[i_penalization_lambda], "is NA"), quote = FALSE)
  } else if (is.null(BICLUST_RESULTS[i_penalization_lambda])) {
    print(paste("penalization_lambda", penalization_lambdas[i_penalization_lambda], "is NULL"), quote = FALSE)
  }else {
    print(paste("penalization_lambda", penalization_lambdas[i_penalization_lambda], "is ok with rand index", BICLUST_RESULTS[[i_penalization_lambda]]$rand_index), quote = FALSE)
  }
}

# print(paste("rand_index for result vs true cluster:", biclust_result$rand_index), quote = FALSE)
# print(paste("Number of iterations:", biclust_result$n_iterations), quote = FALSE)
# print(paste("Silhoutte of first disturbed cluster likelihood (aka how complex was the first likelihood):", biclust_result$db), quote = FALSE)
# print(paste("BIC_all:", biclust_result$BIC), quote = FALSE)
# print(paste("taget_genes_residual_var:"), quote = FALSE)
# print(biclust_result$taget_genes_residual_var, quote = FALSE)