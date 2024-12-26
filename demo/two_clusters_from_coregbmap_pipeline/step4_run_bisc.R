#!/usr/bin/Rscript
rm(list = ls())

library(here)  # To work with paths
library(patchwork)
sink()

# options(warn=2)  # To convert warning messages into error messages which display row of error. For debugging.

# Get absolute path where script is located, by using relative paths.
demo_path <- here::here("demo")
R_path <- here::here("R")
path_data <- here::here('data')
output_path <- demo_path
path_general_env_data <- file.path(path_data, "env_data_general_2clusters_sctransform.RData")
path_cluster1 <- file.path(path_data, "env_data_cluster1.RData")
path_cluster2 <- file.path(path_data, "env_data_cluster2.RData")

source(file.path(R_path, "generate_dummy_data_for_cell_clustering.R"))
source(file.path(R_path, "bisc.R"))
source(file.path(R_path, "randomise_cluster_labels.R"))

# Set seed for example
set.seed(250)
load(path_general_env_data)
load(path_cluster1)
load(path_cluster2)

n_cell_clusters <- 2
n_target_gene_modules <- c(2, 8)  # Number of target gene clusters in each cell cluster
# We assume cells are ordered in the order of cell clusters. So the first x columns are cell cluster 1, etc.
n_cells <- c(ncol(cell_cluster_1), ncol(cell_cluster_2))
true_cluster_allocation <- factor(rep(1:n_cell_clusters, times = n_cells))
rm(n_cells)

penalization_lambdas <- c(0.5)

is_regulator <- inverse_which(indices = ind_reggenes, output_length = n_regulator_genes + n_target_genes)

# for (i in seq(1)) {
#   CELL_CLUSTER_1_RESULTS <- scregclust::scregclust(
#     expression = cell_cluster_1,  # p rows of genes, n columns of cells
#     genesymbols = 1:(n_target_genes + n_regulator_genes),  # Gene row numbers
#     is_regulator = is_regulator,  # Vector indicating which genes are regulators
#     n_modules = n_target_gene_modules[1],
#     penalization = penalization_lambdas,
#     verbose = TRUE,
#     n_cycles = 200,
#     center = FALSE,
#   )
#
#   CELL_CLUSTER_2_RESULTS <- scregclust::scregclust(
#     expression = cell_cluster_2,  # p rows of genes, n columns of cells
#     genesymbols = 1:(n_target_genes + n_regulator_genes),  # Gene row numbers
#     is_regulator = is_regulator,  # Vector indicating which genes are regulators
#     n_modules = n_target_gene_modules[2],
#     penalization = penalization_lambdas,
#     verbose = TRUE,
#     n_cycles = 200,
#     center = FALSE,
#   )
#
#   # Remove genes put in noise cluster by scregclust
#   non_noise_regulator_genes <- vector(length = n_regulator_genes)
#   for (i_target_gene_cluster in seq(2)) {
#     non_noise_regulator_genes <- non_noise_regulator_genes + CELL_CLUSTER_1_RESULTS$results[[1]]$output[[1]]$reg_table[[i_target_gene_cluster]] != 0
#     non_noise_regulator_genes <- non_noise_regulator_genes + CELL_CLUSTER_2_RESULTS$results[[1]]$output[[1]]$reg_table[[i_target_gene_cluster]] != 0
#   }
#   n_regulator_genes <- sum(non_noise_regulator_genes)
#
#   non_noise_genes <- CELL_CLUSTER_1_RESULTS$results[[1]]$output[[1]]$cluster  # Target genes that are in a cluster
#   non_noise_genes <- non_noise_genes + CELL_CLUSTER_2_RESULTS$results[[1]]$output[[1]]$cluster  # Target genes that are in a cluster
#   non_noise_genes <- non_noise_genes != -2  # Convert them to TRUE/FALSE
#   n_target_genes <- sum(non_noise_genes, na.rm = TRUE)
#   ind_targetgenes <- which(c(rep(1, n_target_genes), rep(0, n_regulator_genes)) == 1)
#   ind_reggenes <- which(c(rep(0, n_target_genes), rep(1, n_regulator_genes)) == 1)
#   is_regulator <- inverse_which(indices = ind_reggenes, output_length = n_regulator_genes + n_target_genes)
#
#   non_noise_genes[is.na(non_noise_genes)] <- non_noise_regulator_genes  # add regulator genes that aren't noise
#
#   cell_cluster_1 <- cell_cluster_1[non_noise_genes,]
#   cell_cluster_2 <- cell_cluster_2[non_noise_genes,]
# }


biclust_input_data <- cbind(cell_cluster_1, cell_cluster_2)
initial_clustering <- true_cluster_allocation


all_res <- vector(mode = "list", length = 300)
for (c_seed in seq(300)){
  set.seed(c_seed)
  BICLUST_RESULTS <- vector(mode = "list", length = length(penalization_lambdas))
  sink()
  sink()
  sink()

  for (i_penalization_lambda in seq_along(penalization_lambdas)) {
    print("", quote = FALSE)
    print(paste("Running biclust for penalization_lambda", penalization_lambdas[i_penalization_lambda]), quote = FALSE)
    BICLUST_RESULTS[[i_penalization_lambda]] <- bisc(dat = t(biclust_input_data),
                                                     cell_id = colnames(biclust_input_data),
                                                     true_cell_cluster_allocation = true_cluster_allocation,
                                                     max_iter = 100,
                                                     n_target_gene_clusters = n_target_gene_modules,
                                                     initial_clustering = initial_clustering,
                                                     n_cell_clusters = n_cell_clusters,
                                                     ind_targetgenes = ind_targetgenes,
                                                     ind_reggenes = ind_reggenes,
                                                     output_path = output_path,
                                                     penalization_lambda = penalization_lambdas[i_penalization_lambda],
                                                     calculate_optimization_target = FALSE,
                                                     calculate_silhoutte = FALSE,
                                                     calculate_davies_bouldin_index = FALSE,
                                                     use_garbage_cluster_targets = FALSE)

  }
  all_res[[c_seed]] <- BICLUST_RESULTS

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
