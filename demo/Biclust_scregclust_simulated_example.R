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
path_data <- here::here('data')

source(file.path(R_path, "generate_dummy_data_for_cell_clustering.R"))
source(file.path(R_path, "biclust_scregclust.R"))


#############################################
############ data for dev ###################
#############################################

# Set seed for example
set.seed(1234)


redo_flag = T

# Set variables ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

n_cell_clusters = 2
n_target_gene_clusters = c(10, 7)  # Number of target gene clusters in each cell cluster
n_target_genes = 500             # 2193 from vignette
n_regulator_genes = 100            # 493 from vignette
n_cells = c(1000, 1000)            # c(6000, 6000) from vignette
regulator_means = c(0, 0)         # For generating dummy data, regulator mean in each cell cluster
coefficient_means = list(c(0.0417,
                           0.0343,
                           0.0576,
                           0.043 ,
                           0.0576,
                           0.0413,
                           0.0473,
                           0.0444,
                           0.0481,
                           -0.0139),
                         c(0.0404,
                           0.0519,
                           0.0545,
                           0.0915,
                           0.0663,
                           0.0512,
                           -0.0064
                         )
)  # For generating dummy data, coefficient means in each cell cluster
coefficient_sds = list(c(0.0556,
                         0.037,
                         0.0638,
                         0.0466,
                         0.0761,
                         0.0471,
                         0.0468,
                         0.0611,
                         0.0623,
                         0.0394
),
c(0.0462,
  0.0496,
  0.0807,
  0.1086,
  0.071,
  0.0716,
  0.057
)
)
disturbed_fraction = 0.1  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
plot_stuff = TRUE
plot_suffix = "vignette"
testing_penalization = c(0.1, 0.3) #vector of length n_cell_clusters
# Generate dummy data for each cell cluster that we want ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# generated_data <- generate_dummy_data_for_cell_clustering(
#   n_cell_clusters = n_cell_clusters,
#   n_target_gene_clusters = n_target_gene_clusters,  # Number of target gene clusters in each cell cluster
#   n_target_genes = n_target_genes,
#   n_regulator_genes = n_regulator_genes,
#   n_cells = n_cells,
#   regulator_means = regulator_means,  # For generating dummy data, regulator mean in each cell cluster
#   coefficient_means = coefficient_means,  # For generating dummy data, coefficient means in each cell cluster
#   coefficient_sds = coefficient_sds,
#   disturbed_fraction = 0.25,  # TODO: Add size disturbance too
#   plot_stuff = TRUE,
#  plot_suffix = "vignette"
# )





if (
      !file.exists(file.path(path_data, "env_sim_vgn_data_biclust_sc.rds")) |  redo_flag
    ) {

  generated_data <- generate_dummy_data_for_cell_clustering(
    n_cell_clusters = n_cell_clusters,
    n_target_gene_clusters = n_target_gene_clusters,  # Number of target gene clusters in each cell cluster
    n_target_genes = n_target_genes,          #from vignette
    n_regulator_genes = n_regulator_genes,        # from vignette
    n_cells =n_cells,
    regulator_means = regulator_means,  # For generating dummy data, regulator mean in each cell cluster
    coefficient_means = coefficient_means,
    coefficient_sds = coefficient_sds,
    disturbed_fraction = 0.1,  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
    plot_stuff = TRUE,
    plot_suffix = "vignette",
    testing_penalization = c(0.1, 0.3) #vector of length n_cell_clusters
  )

  saveRDS(generated_data, file.path(path_data, "env_sim_vgn_data_biclust_sc.rds"))

} else{

  generated_data <- readRDS(file.path(path_data, "env_sim_vgn_data_biclust_sc.rds"))

}



# Because "dat <- cbind(Z_t, Z_r)" in generate_dummy_data_for_cell_clustering
ind_targetgenes <- which(c(rep(1, n_target_genes), rep(0, n_regulator_genes)) == 1)
ind_reggenes <- which(c(rep(0, n_target_genes), rep(1, n_regulator_genes)) == 1)


disturbed_initial_cell_clust <- factor(generated_data$disturbed_initial_cell_clust)

biclust_input_data <- generated_data$dat
colnames(biclust_input_data) <- c(paste0("t", 1:n_target_genes), paste0("r", 1:n_regulator_genes))
biclust_input_data <- tibble::as_tibble(biclust_input_data)

# # These needs to be strings for discrete labels in pca plot
# data_for_plotting <- tibble::as_tibble(true_cell_ uster_allocation = generated_data$true_cell_clust,
#                                        biclust_input_data)
# pca_res <- prcomp(biclust_input_data, scale. = TRUE)}
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

# Setup variables that will be used throughout the script
# We assume target genes start with t then a number. And likewise for regulator genes.
n_total_cells_train <- length(train_indices)
n_total_cells <- nrow(biclust_input_data)
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


penalization_lambdas <- c( 0.5) # c( 0.00001, 0.1, 0.2, 0.5)
BICLUST_RESULTS <- vector(mode = "list", length = length(penalization_lambdas))

max_iter <- 20

if (
      !file.exists(file.path(path_data, "env_sim_vgn_res_biclust_sc.rds"))  |  redo_flag
    ) {

  for (i_penalization_lambda in seq_along(penalization_lambdas)) {
    print("", quote = FALSE)
    print(paste("Running biclust for penalization_lambda", penalization_lambdas[i_penalization_lambda]), quote = FALSE)
    BICLUST_RESULTS[[i_penalization_lambda]] <- biclust_scregclust(
      dat = biclust_input_data,
      cell_id = cell_id,
      true_cell_cluster_allocation = factor(generated_data$true_cell_clust),
      max_iter = max_iter,
      n_target_gene_clusters,
      initial_clustering,
      n_cell_clusters,
      ind_targetgenes,
      ind_reggenes,
      output_path,
      penalization_lambda = penalization_lambdas[i_penalization_lambda],
      use_complex_cluster_allocation = FALSE,
      calculate_BIC = FALSE,
      calculate_silhoutte = FALSE,
      calculate_davies_bouldin_index = FALSE,
      plot_suffix = "vignette",
      always_use_flat_prior = FALSE,
      use_garbage_cluster_targets  = F
    )
  }

  saveRDS(BICLUST_RESULTS, file.path(path_data, "env_sim_vgn_res_biclust_sc.rds"))

} else{

  BICLUST_RESULTS <- readRDS(file.path(path_data, "env_sim_vgn_res_biclust_sc.rds"))

}

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
