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


#############################################
############ data for dev ###################
#############################################

# Set seed for example
set.seed(214)

# Set variables ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

n_cell_clusters <- 2
n_target_gene_clusters <- c(2, 3)  # Number of target gene clusters in each cell cluster
n_target_genes <- 30
n_regulator_genes <- 70
n_cells <- c(10000, 10000)
regulator_means <- c(5, 1)  # For generating dummy data, regulator mean in each cell cluster
coefficient_means <- list(c(1, 20), c(5, 20, 100))  # For generating dummy data, coefficient means in each cell cluster
coefficient_sds <- list(c(0.1, 0.1), c(0.1, 0.1, 0.1))
true_cluster_allocation <- rep(1:n_cell_clusters, times = n_cells)
n_total_cells <- sum(n_cells)

# Generate dummy data for each cell cluster that we want ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
generated_data <- generate_dummy_data_for_cell_clustering(
  n_cell_clusters = n_cell_clusters,
  n_target_gene_clusters = n_target_gene_clusters,  # Number of target gene clusters in each cell cluster
  n_target_genes = n_target_genes,
  n_regulator_genes = n_regulator_genes,
  n_cells = n_cells,
  regulator_means = regulator_means,  # For generating dummy data, regulator mean in each cell cluster
  coefficient_means = coefficient_means,  # For generating dummy data, coefficient means in each cell cluster
  coefficient_sds = coefficient_sds,
  disturbed_fraction = 0.49  # TODO: Add size disturbance too
)


# Because "dat <- cbind(Z_t, Z_r)" in generate_dummy_data_for_cell_clustering
ind_targetgenes <- which(c(rep(1, n_target_genes), rep(0, n_regulator_genes)) == 1)
ind_reggenes <- which(c(rep(0, n_target_genes), rep(1, n_regulator_genes)) == 1)


disturbed_initial_cell_clust <- factor(generated_data$disturbed_initial_cell_clust)

biclust_input_data <- generated_data$dat
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

max_iter <- 50
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

penalization_lambdas <- c(0.0001, 0.001, 0.01, 0.1, 0.2, 0.3, 0.5, 0.7)
BICLUST_RESULTS <- vector(mode = "list", length = length(penalization_lambdas))


for (i_penalization_lambda in seq_along(penalization_lambdas)) {
  print("", quote = FALSE)
  print(paste("Running biclust for penalization_lambda", penalization_lambdas[i_penalization_lambda]), quote = FALSE)
  BICLUST_RESULTS[[i_penalization_lambda]] <- biclust(dat = biclust_input_data,
                                                      cell_id = cell_id,
                                                      true_cell_cluster_allocation = factor(generated_data$true_cell_clust),
                                                      max_iter = 50,
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