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
set.seed(2147)

# Load data
# Folders and filenames
path_data <- here::here('data')
path_medium <- file.path(path_data, "medium.rds")
path_medium_cleaned <- file.path(path_data, "medium_cleaned.rds")
path_cell_types <- file.path(path_data, "cell_types.rds")
path_env_data <- file.path(path_data, "env_data.RData")
# Data from https://cellxgene.cziscience.com/collections/999f2a15-3d7e-440b-96ae-2c806799c08c
if (file.exists(path_env_data)) {
  load(path_env_data)
}else {
  if (file.exists(path_medium_cleaned) && file.exists(path_cell_types)) {
    # Read data
    d <- readRDS(path_medium_cleaned)
    cell_types <- readRDS(path_cell_types)
  }else {
    # Read data
    d <- readRDS(path_medium)

    # Keep the correct type of cells
    cell_types <- d@meta.data$annotation_level_3
    keep_cells <- cell_types == 'AC-like' |
      cell_types == 'MES-like' |
      cell_types == 'NPC-like' |
      cell_types == 'OPC-like'

    keep_cells_split <- sample(c(1, 2), length(keep_cells), prob = c(0.25, 0.75), replace = T)
    keep_ind <- keep_cells_split == 1
    keep_cells <- keep_cells & keep_ind

    cell_types <- factor(cell_types[keep_cells])
    saveRDS(cell_types, path_cell_types)
    d <- d[, keep_cells]

    # Remove genes which are 'constant', aka no variance, just one number.
    # scregclust doesn't work on those.
    d_temp <- d[["RNA"]]$data
    non_constant_ind <- which(apply(d_temp, 1, sd) != 0)
    d <- d[non_constant_ind,]
    rm(non_constant_ind)
    rm(d_temp)
    saveRDS(d, path_medium_cleaned)
  }

  # Find out which are regulator genes
  # Make it into a matrix and rename row names from
  # Ensembl gene IDs (like "ENSG00000269696") to standard gene names (like "A1BG"),
  fake_matrix <- matrix(0, nrow = nrow(d), ncol = 1)  # This is fast
  # dm <- tibble::as_tibble(as.matrix(GetAssayData(d, assay = "RNA", slot = "counts")))  # Using the data is slow
  rownames(fake_matrix) <- as.vector(d@assays$RNA@meta.features$feature_name)
  out <- scregclust::scregclust_format(fake_matrix)  # Needs to be a matrix to use this
  is_regulator <- out[['is_regulator']]


  # Put some stuff in order gene/row wise.
  # This puts the regulator genes at the end
  d <- d[order(is_regulator),]
  is_regulator <- is_regulator[order(is_regulator)]

  # Sub sample genes
  n_target_genes <- 300
  n_regulator_genes <- 150
  old_ind_regulators <- which(as.logical(is_regulator))
  new_ind_targetgenes <- sample(1:(min(old_ind_regulators) - 1), n_target_genes)
  new_ind_regulators <- sample(old_ind_regulators, n_regulator_genes)
  new_ind_genes <- sort(c(new_ind_targetgenes, new_ind_regulators))
  keep_genes <- replace(rep(FALSE, length(is_regulator)), new_ind_genes, TRUE)
  #  Apply sub sampling
  d <- d[keep_genes,]
  is_regulator <- (c(rep(0, n_target_genes), rep(1, n_regulator_genes)) == 1)

  # Put stuff in order cell/column wise
  cell_order <- order(cell_types)
  cell_types <- cell_types[cell_order]
  unique(cell_types[cell_order])
  d <- d[, cell_order]

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
                                                                  fraction_randomised = 0.10))
  # convert from seurat object into matrix that we can work with
  # List of object layers
  # SeuratObject::Layers(d)
  biclust_input_data <- as.matrix(Seurat::GetAssayData(d, assay = "RNA", layer = "data"))
  rm(d)
  biclust_input_data <- t(tibble::as_tibble(biclust_input_data))
  colnames(biclust_input_data) <- c(paste0("t", 1:n_target_genes), paste0("r", 1:n_regulator_genes))


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

  max_iter <- 10
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

  penalization_lambdas <- c(0.001)
  BICLUST_RESULTS <- vector(mode = "list", length = length(penalization_lambdas))
}



###########################################
####initialise variables for dev ##########
##########################################

max_iter <- 100
initial_clustering <- disturbed_initial_cell_clust
# n_target_genes <- n_target_genes
# n_regulator_genes <- n_regulator_genes
# n_total_cells <- n_total_cells
# n_cell_clusters <- n_cell_clusters
# ind_targetgenes <- ind_targetgenes
# ind_reggenes <- ind_reggenes

# reg_split <- sample(c(1, 2), length(ind_reggenes), prob = c(0.8, 0.2), replace = T)
# remove_indices <- which(cell_data_split == 1)
# remove_ind_reggenes <- ind_reggenes[remove_indices]
# ind_reggenes <- ind_reggenes[-remove_ind_reggenes]
#

# output_path <- modded_output_path

i_cell_cluster <- 1

use_weights <- FALSE
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


penalization_lambda <- 0.001
results_outer <- vector("list", 4)

# cell_cluster <- 1
# iter <- 1

library(scregclust)
for(cell_cluster in seq(1,4)){
  print(paste("calculating for cell cluster", cell_cluster))
  current_cell_cluster <- biclust_input_data[true_cluster_allocation == cell_cluster,]

  #run screg with a bunch of different cluster setings
  results_inner <- vector("list", 7)

  for(iter in seq(1,7)){
    target_gene_clusters <- seq(2,8)[iter]

    print(paste("calculating for",target_gene_clusters, "clusters"))


    current_cell_cluster <- current_cell_cluster[,
                                             apply(current_cell_cluster,
                                                     MARGIN = 2,
                                                   function(x) sum(x != 0)) > 200]


    n_target_genes <- length(grep("^t", colnames(current_cell_cluster), value = TRUE))
    n_regulator_genes <- length(grep("^r", colnames(current_cell_cluster), value = TRUE))
    ind_targetgenes <- which(c(rep(1, n_target_genes), rep(0, n_regulator_genes)) == 1)
    ind_reggenes <- which(c(rep(0, n_target_genes), rep(1, n_regulator_genes)) == 1)
    colnames(current_cell_cluster) <- c(paste0("t", 1:n_target_genes), paste0("r", 1:n_regulator_genes))




    results_inner[[iter]] <- scregclust::scregclust(
      expression = t(current_cell_cluster),  # p rows of genes, n columns of cells
      split_indices = NULL,
      genesymbols = 1:(n_target_genes + n_regulator_genes),  # Gene row numbers
      is_regulator = inverse_which(indices = ind_reggenes, output_length = n_regulator_genes + n_target_genes),  # Vector indicating which genes are regulators
      n_cl = target_gene_clusters,
      penalization = penalization_lambda,
      noise_threshold = 0.00001,
      verbose = TRUE,
      n_cycles = 150,
      compute_silhouette = T
    )
  }

  results_outer[[cell_cluster]] <- results_inner

  saveRDS(results_inner, file.path(path_data,paste0("screg_results_cluster_",cell_cluster,
                                              "_penalization_",penalization_lambda,
                                              ".rds")))
}


for(cell_cluster in seq(1,4)){
  print(paste("calculating for cell cluster", cell_cluster))

  #run screg with a bunch of different cluster setings
  results_inner <- vector("list", 7)

  for(iter in seq(1,7)){
    target_gene_clusters <- seq(2,8)[iter]
    print(file.path(path_data,paste0("screg_results_cluster_",cell_cluster,
                                     "_penalization_",penalization_lambda,
                                     ".rds"))
          )
  }
}


# scregclust::plot_silhouettes(list_of_fits = results_outer[[1]],
                             # penalization = c(penalization_lambda))

demo_path <- here::here("demo")
R_path <- here::here("R")
output_path <- demo_path

for(cell_cluster in seq(1,4)) {

  png(
  file.path(output_path,paste0("cluster_count_helper_coreGB_cell_cluster",
                               cell_cluster,
                               "_penalisation_",
                               penalization_lambda,
                               ".png")
            )
  ,width = 1024, height = 480, units = "px"
      )

 scregclust::plot_cluster_count_helper(
          list_of_fits = results_outer[[cell_cluster]],
          penalization = c(penalization_lambda)
          ) -> p
  print(p)
  dev.off()

}


