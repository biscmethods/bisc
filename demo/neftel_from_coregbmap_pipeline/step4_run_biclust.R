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
path_coreGBmap <- file.path(path_data, "medium.rds")
path_env_data_neftel2019 <- file.path(path_data, "env_data_findtargetcluster_sctransform_neftel2019.RData")
path_general_env_data_neftel2019 <- file.path(path_data, "env_data_general_findtargetcluster_sctransform_neftel2019.RData")
path_AC_neftel2019 <- file.path(path_data, "env_data_AC_sctransform_neftel2019.RData")
path_MES_neftel2019 <- file.path(path_data, "env_data_MES_sctransform_neftel2019.RData")
path_NPC_neftel2019 <- file.path(path_data, "env_data_NPC_sctransform_neftel2019.RData")
path_OPC_neftel2019 <- file.path(path_data, "env_data_OPC_sctransform_neftel2019.RData")

source(file.path(R_path, "generate_dummy_data_for_cell_clustering.R"))
source(file.path(R_path, "biclust_scregclust.R"))
source(file.path(R_path, "randomise_cluster_labels.R"))

# Set seed for example
set.seed(250)
load(path_general_env_data_neftel2019)
load(path_AC_neftel2019)
load(path_MES_neftel2019)

n_cell_clusters <- 2
n_target_gene_clusters <- c(2, 2)  # Number of target gene clusters in each cell cluster
# We assume cells are ordered in the order of cell clusters. So the first x columns are cell cluster 1, etc.
n_cells <- c(ncol(cell_cluster_AC), ncol(cell_cluster_MES))
true_cluster_allocation <- factor(rep(1:n_cell_clusters, times = n_cells))
rm(n_cells)

penalization_lambdas <- c(0.2)

is_regulator <- inverse_which(indices = ind_reggenes, output_length = n_regulator_genes + n_target_genes)

for (i in seq(1)) {
  AC <- scregclust::scregclust(
    expression = cell_cluster_AC,  # p rows of genes, n columns of cells
    genesymbols = 1:(n_target_genes + n_regulator_genes),  # Gene row numbers
    is_regulator = is_regulator,  # Vector indicating which genes are regulators
    n_cl = 2,
    penalization = penalization_lambdas,
    verbose = TRUE,
    n_cycles = 200,
    center=FALSE,
  )

  MES <- scregclust::scregclust(
    expression = cell_cluster_MES,  # p rows of genes, n columns of cells
    genesymbols = 1:(n_target_genes + n_regulator_genes),  # Gene row numbers
    is_regulator = is_regulator,  # Vector indicating which genes are regulators
    n_cl = 2,
    penalization = penalization_lambdas,
    verbose = TRUE,
    n_cycles = 200,
    center=FALSE,
  )

  # Remove genes put in noise cluster by scregclust
  non_noise_regulator_genes <- vector(length = n_regulator_genes)
  for (i_target_gene_cluster in seq(2)) {
    non_noise_regulator_genes <- non_noise_regulator_genes + AC$results[[1]]$output[[1]]$reg_table[[i_target_gene_cluster]] != 0
    non_noise_regulator_genes <- non_noise_regulator_genes + MES$results[[1]]$output[[1]]$reg_table[[i_target_gene_cluster]] != 0
  }
  n_regulator_genes <- sum(non_noise_regulator_genes)

  non_noise_genes <- AC$results[[1]]$output[[1]]$cluster  # Target genes that are in a cluster
  non_noise_genes <- non_noise_genes + MES$results[[1]]$output[[1]]$cluster  # Target genes that are in a cluster
  non_noise_genes <- non_noise_genes != -2  # Convert them to TRUE/FALSE
  n_target_genes <- sum(non_noise_genes, na.rm = TRUE)
  ind_targetgenes <- which(c(rep(1, n_target_genes), rep(0, n_regulator_genes)) == 1)
  ind_reggenes <- which(c(rep(0, n_target_genes), rep(1, n_regulator_genes)) == 1)
  is_regulator <- inverse_which(indices = ind_reggenes, output_length = n_regulator_genes + n_target_genes)

  non_noise_genes[is.na(non_noise_genes)] <- non_noise_regulator_genes  # add regulator genes that aren't noise

  cell_cluster_AC <- cell_cluster_AC[non_noise_genes,]
  cell_cluster_MES <- cell_cluster_MES[non_noise_genes,]
}


biclust_input_data <- cbind(cell_cluster_AC, cell_cluster_MES)
initial_clustering <- (true_cluster_allocation)


all_res <- vector(mode = "list", length = 300)
for (c_seed in seq(300)){
  set.seed(c_seed)
  BICLUST_RESULTS <- vector(mode = "list", length = length(penalization_lambdas))

  for (i_penalization_lambda in seq_along(penalization_lambdas)) {
    print("", quote = FALSE)
    print(paste("Running biclust for penalization_lambda", penalization_lambdas[i_penalization_lambda]), quote = FALSE)
    BICLUST_RESULTS[[i_penalization_lambda]] <- biclust(dat = t(biclust_input_data),
                                                        cell_id = colnames(biclust_input_data),
                                                        true_cell_cluster_allocation = true_cluster_allocation,
                                                        max_iter = 100,
                                                        n_target_gene_clusters = n_target_gene_clusters,
                                                        initial_clustering = initial_clustering,
                                                        n_cell_clusters = n_cell_clusters,
                                                        ind_targetgenes = ind_targetgenes,
                                                        ind_reggenes = ind_reggenes,
                                                        output_path = output_path,
                                                        penalization_lambda = penalization_lambdas[i_penalization_lambda],
                                                        use_complex_cluster_allocation = TRUE,
                                                        calculate_BIC = FALSE,
                                                        calculate_silhoutte = FALSE,
                                                        calculate_davies_bouldin_index = FALSE)

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
