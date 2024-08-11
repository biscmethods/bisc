  #!/usr/bin/Rscript

library(tidyverse)
library(aricode)     # To calculate rand index
library(ggplot2)     # To plot things #TODO: What things? literally anything
library(ggalluvial)  # To plot thingsv #TODO: What things?
library(reshape2)
library(here)        # To work with paths
library(ggfortify)   # For pca-plot
library(pracma)      # For pseudo inverse
library(stats)       # For kmeans
library(ppclust)
library(clusterSim)  # for db

# Get absolute path where script is located, by using relative paths.
R_path <- here::here("R")
source(file.path(R_path, "plot_cluster_history.R"))
source(file.path(R_path, "plot_loglikelihood.R"))
source(file.path(R_path, "generate_dummy_data_for_cell_clustering.R"))

  # Set seed for example
  set.seed(231)

  # Set variables ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # n_cell_clusters <- 2
  # n_target_gene_clusters <- c(2, 3)  # Number of target gene clusters in each cell cluster
  # n_target_genes <- 50
  # n_regulator_genes <- 30
  # n_cells <- c(10000, 5000)
  # regulator_means <- c(5, 1)  # For generating dummy data, regulator mean in each cell cluster
  # coefficient_means <- list(c(1, 20), c(5, 20, 100))  # For generating dummy data, coefficient means in each cell cluster
  # coefficient_sds <- list(c(0.1, 0.1), c(0.1, 0.1, 0.1))
  # true_cluster_allocation <- rep(1:n_cell_clusters, times = n_cells)
  # n_total_cells <- sum(n_cells)

  # Generate dummy data for each cell cluster that we want ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # generated_data <- generate_dummy_data_for_cell_clustering(
    # n_cell_clusters = n_cell_clusters,
    # n_target_gene_clusters = n_target_gene_clusters,  # Number of target gene clusters in each cell cluster
    # n_target_genes = n_target_genes,
    # n_regulator_genes = n_regulator_genes,
    # n_cells = n_cells,
    # regulator_means = regulator_means,  # For generating dummy data, regulator mean in each cell cluster
    # coefficient_means = coefficient_means,  # For generating dummy data, coefficient means in each cell cluster
    # coefficient_sds = coefficient_sds,
    # disturbed_fraction = 0.22
  # )

  generated_data <- generate_dummy_data_for_cell_clustering(
  )

  str(generated_data)

  ind_reggenes <- which(c(rep(0, n_target_genes), rep(1, n_regulator_genes)) == 1)
  ind_targetgenes <- which(c(rep(1, n_target_genes), rep(0, n_regulator_genes)) == 1)

  disturbed_initial_cell_clust <- factor(generated_data$disturbed_initial_cell_clust)

  biclust_input_data <- generated_data$dat
  colnames(biclust_input_data) <- c(paste0("t", 1:n_target_genes), paste0("r", 1:n_regulator_genes))
  biclust_input_data <- tibble::as_tibble(biclust_input_data)

  example_regulator <- biclust_input_data[,ind_reggenes[1]]

  example_regulator %>% ggplot()

  betas <- generated_data$true_betas
  str(generated_data$dat)
  str(generated_data$true_cell_clust)

  hist(betas[[1]][[1]])
  hist(betas[[1]][[2]])

  hist(betas[[2]][[1]])
  hist(betas[[2]][[2]])
  hist(betas[[2]][[3]])

  #ta en regulatorgen och en targetgen och se hur vägen till targetgen ser ut
