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

path_env_data <- file.path(path_data, "env_data_pbmc_sctransform.RData")


source(file.path(R_path, "generate_dummy_data_for_cell_clustering.R"))
source(file.path(R_path, "bisc.R"))
source(file.path(R_path, "randomise_cluster_labels.R"))

# Set seed for example
set.seed(250)
load(path_env_data)

keep_cluster <- c(1, 2, 3, 5, 7)
n_target_gene_modules <- c(2, 3, 3, 5, 2)  # Number of target gene clusters in each cell cluster
n_cell_clusters <- length(keep_cluster)

logical_keep_cell_cluster <- true_cluster_allocation == 100
for(i_cell_cluster in keep_cluster){
  logical_keep_cell_cluster <- logical_keep_cell_cluster | (true_cluster_allocation==i_cell_cluster)
}

biclust_input_data <- d[,which(logical_keep_cell_cluster)]
true_cluster_allocation <- true_cluster_allocation[which(logical_keep_cell_cluster)]
true_cluster_allocation <- match(true_cluster_allocation, sort(unique(true_cluster_allocation)))


penalization_lambdas <- c(0.5)
initial_clustering <- factor(sample(true_cluster_allocation))


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
                                                     max_iter = 10,
                                                     n_target_gene_clusters = n_target_gene_modules,
                                                     initial_clustering = initial_clustering,
                                                     n_cell_clusters = n_cell_clusters,
                                                     ind_targetgenes = ind_targetgenes,
                                                     ind_reggenes = ind_reggenes,
                                                     output_path = output_path,
                                                     penalization_lambda = penalization_lambdas[i_penalization_lambda],
                                                     use_complex_cluster_allocation = FALSE,
                                                     calculate_BIC = FALSE,
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
