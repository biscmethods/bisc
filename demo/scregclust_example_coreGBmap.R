#!/usr/bin/Rscript
rm(list = ls())

library(here)  # To work with paths
library(patchwork)
sink()
#TODO: Subsample genes?
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
path_medium <- file.path(path_data, "medium.rds")  # Data from https://cellxgene.cziscience.com/collections/999f2a15-3d7e-440b-96ae-2c806799c08c

# Read data
d <- readRDS(path_medium)


# Keep the correct type of cells
cell_types <- d@meta.data$annotation_level_3
keep_cells <- cell_types == 'AC-like' |
  cell_types == 'MES-like' |
  cell_types == 'NPC-like' |
  cell_types == 'OPC-like'

keep_cells_split <- sample(c(1, 2), length(keep_cells), prob = c(0.8, 0.2), replace = T)
keep_ind <- keep_cells_split == 1
keep_cells <- keep_cells & keep_ind

cell_types <- factor(cell_types[keep_cells])
d <- d[,keep_cells]

# Find out which are regulator genes
# Make it into a matrix and rename row names from
# Ensembl gene IDs (like "ENSG00000269696") to standard gene names (like "A1BG"),
fake_matrix <- matrix(0, nrow = nrow(d), ncol = 1)  # This is fast
# dm <- tibble::as_tibble(as.matrix(GetAssayData(d, assay = "RNA", slot = "counts")))  # Using the data is slow
rownames(fake_matrix) <- as.vector(d@assays$RNA@meta.features$feature_name)
out <- scregclust::scregclust_format(fake_matrix)  # Needs to be a matrix to use this
is_regulator <- out[['is_regulator']]

# Put some stuff in order gene/row wise
d <- d[order(is_regulator),]
is_regulator <- is_regulator[order(is_regulator)]

# Put stuff in order cell/column wise
cell_order <- order(cell_types)
cell_types <- cell_types[cell_order]
unique(cell_types[cell_order])
d <- d[,cell_order]

# Now when we have put the cells in order we just need to count the cells
# in each cell cluster. Then it's easy to make a vector with the true_cluster_allocation
# further down in code
n_cells_cell_cluster_1 <- sum(cell_types == 'AC-like')
n_cells_cell_cluster_2 <- sum(cell_types == 'MES-like')
n_cells_cell_cluster_3 <- sum(cell_types == 'NPC-like')
n_cells_cell_cluster_4 <- sum(cell_types == 'OPC-like')

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
total_n_cells_in_genedata <- ncol(d)
print(paste(str(total_n_cells_in_cellstxt), " = ", str(total_n_cells_in_genedata)), quote = FALSE)


# Set variables ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

n_cell_clusters <- 4
n_target_gene_clusters <- c(3, 3, 3, 3)  # TODO:Number of target gene clusters in each cell cluster
n_target_genes <- length(is_regulator) - sum(is_regulator == 1)
n_regulator_genes <- sum(is_regulator == 1)
n_total_cells <- ncol(d)
# We assume cells are ordered in the order of cell clusters. So the first x columns are cell cluster 1, etc.
n_cells <- c(n_cells_cell_cluster_1, n_cells_cell_cluster_2, n_cells_cell_cluster_3, n_cells_cell_cluster_4)
true_cluster_allocation <- rep(1:n_cell_clusters, times = n_cells)  # TODO: do this


# Because "dat <- cbind(Z_t, Z_r)" in generate_dummy_data_for_cell_clustering
ind_targetgenes <- which(c(rep(1, n_target_genes), rep(0, n_regulator_genes)) == 1)
ind_reggenes <- which(c(rep(0, n_target_genes), rep(1, n_regulator_genes)) == 1)

disturbed_initial_cell_clust <- factor(randomise_cluster_labels(cluster_labels = true_cluster_allocation,
                                                                fraction_randomised = 0.25))

biclust_input_data <- t(d)
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

# reg_split <- sample(c(1, 2), length(ind_reggenes), prob = c(0.8, 0.2), replace = T)
# remove_indices <- which(cell_data_split == 1)
# remove_ind_reggenes <- ind_reggenes[remove_indices]
# ind_reggenes <- ind_reggenes[-remove_ind_reggenes]
#

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