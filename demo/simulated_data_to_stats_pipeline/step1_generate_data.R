# !/usr/bin/Rscript

library(patchwork)
library(rasterVis)
library(cowplot)
library(grid)
sink()

# options(warn=2)  # To convert warning messages into error messages which display row of error. For debugging.

source(file.path(R_path, "generate_dummy_data_for_cell_clustering.R"))


construct_biclust_matrix <- function(cell_cluster_allocation,
                                     target_gene_allocation){
  n_total_cells <- length(cell_cluster_allocation)
  n_target_genes <- length(target_gene_allocation[[1]])
  # Construct correct biclust matrix --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  true_biclust_matrix <- matrix(0, nrow = n_total_cells, ncol = n_target_genes)

  # Check each entry in the matrix (each cell and gene pair), and assign a string-number to each unique "cell-cluster-gene-module".
  for (i in 1:n_total_cells) {
    cluster <- cell_cluster_allocation[i]
    gene_allocation <- target_gene_allocation[[cluster]][1:n_target_genes]
    gene_allocation[gene_allocation==-1] <- 0
    # Create unique numbers for each pair
    true_biclust_matrix[i, ] <- paste0(cluster, gene_allocation)
  }

  # Convert the string-numbers to numeric matrix (starting from 1 this time)
  true_biclust_matrix <- matrix(as.numeric(true_biclust_matrix), nrow = n_total_cells, ncol = n_target_genes)
  true_biclust_matrix <- matrix(as.integer(as.factor(true_biclust_matrix)), nrow=nrow(true_biclust_matrix), ncol=ncol(true_biclust_matrix))
  # correct_clustering <- as.vector(true_biclust_matrix)
  return(true_biclust_matrix)
}

split_traintest <- function(generated_data,#  generated_data               :List of 9
                            # ..$ disturbed_initial_cell_clust: int [1:400] 1 1 1 1 1 1 1 1 1 1 ...
                            # ..$ initial_cell_clust          : int [1:400] 1 1 1 1 1 1 1 1 1 1 ...
                            # ..$ true_cell_clust             : int [1:400] 1 1 1 1 1 1 1 1 1 1 ...
                            # ..$ true_betas                  :List of 2
                            # .. ..$ :List of 4
                            # .. .. ..$ : num [1:6, 1:100] 0.994 0 0 0 0 ...
                            # .. .. ..$ : num [1:6, 1:100] 0 0 0 0 0 0 0 0 0 0 ...
                            # .. .. ..$ : num [1:6, 1:100] 0 0 0 0 0 0 0 0 0 0 ...
                            # .. .. ..$ : num [1:6, 1:100] 0 0 0 0 0 0 0 0 0 0 ...
                            # .. ..$ :List of 2
                            # .. .. ..$ : num [1:6, 1:100] 0 0 0 0 9.99 ...
                            # .. .. ..$ : num [1:6, 1:100] 0 0 0 0 0 0 0 0 0 0 ...
                            # ..$ dat                         : num [1:400, 1:106] -2.19 -3.14 1.58 -2.82 -2.51 ...
                            # ..$ true_target_gene_allocation :List of 2
                            # .. ..$ : int [1:100] 1 1 1 1 1 1 1 1 1 1 ...
                            # .. ..$ : int [1:100] 1 1 1 1 1 1 1 1 1 1 ...
                            # ..$ true_Pi                     :List of 2
                            # .. ..$ : num [1:4, 1:100] 1 0 0 0 1 0 0 0 1 0 ...
                            # .. ..$ : num [1:2, 1:100] 1 0 1 0 1 0 1 0 1 0 ...
                            # ..$ true_S                      :List of 2
                            # .. ..$ : num [1:4, 1:6] 1 0 0 0 0 1 0 0 0 0 ...
                            # .. ..$ : num [1:2, 1:6] 0 0 0 0 0 0 0 0 1 0 ...
                            # ..$ counts                      : logi NA
                            disturbed_initial_cell_clust,
                            biclust_input_data, # tibble [400 × 106] (S3: tbl_df/tbl/data.frame)
                            n_cells, #  : num [1:2] 200 200
                            n_total_cells,  # num 400
                            cell_id,  # num [1:400] 1 2 3 4 5 6 7 8 9 10 ...
                            correct_clustering
){
  # Split data into train/test
  cell_data_split <- sample(c(1, 2), nrow(biclust_input_data), prob = c(0.5, 0.5), replace = T)
  train_indices <- which(cell_data_split == 1)
  test_indices <- which(cell_data_split == 2)

  biclust_input_data_train <- biclust_input_data[train_indices,]
  biclust_input_data_test <- biclust_input_data[test_indices,]

  # Setup variables that will be used throughout the script
  # We assume target genes start with t then a number. And likewise for regulator genes.
  n_total_cells_train <- length(train_indices)
  n_total_cells_test <- length(test_indices)
  cell_id_train <- cell_id[train_indices]
  cell_id_test <- cell_id[test_indices]

  disturbed_initial_cell_clust_train <- disturbed_initial_cell_clust[train_indices]
  disturbed_initial_cell_clust_test <- disturbed_initial_cell_clust[test_indices]

  initial_clustering_train <- disturbed_initial_cell_clust_train
  initial_clustering_test <- disturbed_initial_cell_clust_test

  # Set up some variables

  true_cell_cluster_allication_train <- generated_data$true_cell_clust[train_indices]
  true_cell_cluster_allication_test <- generated_data$true_cell_clust[test_indices]

  n_cells_train <- as.vector(table(true_cell_cluster_allication_train))
  n_cells_test <- as.vector(table(true_cell_cluster_allication_test))

  # Set output of training data -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  generated_data$data <- as.matrix(biclust_input_data_train)
  generated_data$true_cell_clust <- true_cell_cluster_allication_train
  generated_data$initial_cell_clust <- initial_clustering_train
  generated_data$disturbed_initial_cell_clust <- disturbed_initial_cell_clust_train
  disturbed_initial_cell_clust <- disturbed_initial_cell_clust_train
  biclust_input_data <- biclust_input_data_train
  n_total_cells <- n_total_cells_train
  n_cells <- n_cells_train
  cell_id <- cell_id_train
  correct_clustering <- as.vector(construct_biclust_matrix(cell_cluster_allocation = generated_data$true_cell_clust,
                                                           target_gene_allocation= generated_data$true_target_gene_allocation))


  # Output of test data ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # generated_data$data <- as.matrix(biclust_input_data_train)
  generated_data$true_cell_clust_test <- true_cell_cluster_allication_test
  generated_data$initial_cell_clust_test <- initial_clustering_test
  generated_data$disturbed_initial_cell_clust_test <- disturbed_initial_cell_clust_test
  disturbed_initial_cell_clust_test <- disturbed_initial_cell_clust_test
  biclust_input_data_test <- biclust_input_data_test
  n_total_cells_test <- n_total_cells_test
  n_cells_test <- n_cells_test
  cell_id_test <- cell_id_test
  correct_clustering_test <- as.vector(construct_biclust_matrix(cell_cluster_allocation = generated_data$true_cell_clust_test,
                                                                target_gene_allocation= generated_data$true_target_gene_allocation))

  return(list(generated_data = generated_data,
              disturbed_initial_cell_clust = disturbed_initial_cell_clust,
              disturbed_initial_cell_clust_test = disturbed_initial_cell_clust_test,
              biclust_input_data = biclust_input_data,
              biclust_input_data_test = biclust_input_data_test,
              n_total_cells = n_total_cells,
              n_total_cells_test = n_total_cells_test,
              cell_id = cell_id,
              cell_id_test = cell_id_test,
              correct_clustering = correct_clustering,
              correct_clustering_test = correct_clustering_test,
              n_cells = n_cells,
              n_cells_test = n_cells_test
  ))



}

# Function that creates data for a scenario and returns a list of all variables needed for both generating the same data and run analysis -------------------------------------------------------------------------------------------------
create_scenario <- function(scenario_list,
                            description,
                            plot_heatmap = FALSE,
                            plot_suffix = "",
                            n_cell_clusters,
                            n_target_gene_clusters,
                            n_target_genes,
                            n_regulator_genes,
                            n_cells,
                            regulator_means,
                            coefficient_means,
                            coefficient_sds,
                            disturbed_fraction,
                            testing_penalization_data_gen,
                            trivial_regulator_networks,
                            seed){
  # Set seed for example
  set.seed(seed)


  # This function adds all data that is needed to both use and create a scenario
  add_scenario_data <- function(scenario_list,
                                description,
                                n_cell_clusters,
                                n_target_gene_clusters,
                                n_target_genes,
                                n_regulator_genes,
                                n_cells,
                                n_cells_test,
                                regulator_means,
                                coefficient_means,
                                coefficient_sds,
                                disturbed_fraction,
                                testing_penalization_data_gen,
                                generated_data,
                                ind_targetgenes,
                                ind_reggenes,
                                disturbed_initial_cell_clust,
                                disturbed_initial_cell_clust_test,
                                biclust_input_data,
                                biclust_input_data_test,
                                n_total_cells,
                                n_total_cells_test,
                                cell_id,
                                cell_id_test,
                                correct_clustering,
                                correct_clustering_test,
                                trivial_regulator_networks,
                                seed){
    # Check that all arguments have been given
    args <- names(formals())
    missing_args <- sapply(args, function(x) {
      is.symbol(formals()[[x]]) && missing(as.name(x))
    })
    if(any(missing_args)) {
      stop(paste("Missing arguments:", paste(names(missing_args)[missing_args], collapse=", ")))
    }

    # Create the new scenario entry
    scenario_list <- append(scenario_list,
                            list(list(description = description,
                                      n_cell_clusters = n_cell_clusters,
                                      n_target_gene_clusters = n_target_gene_clusters,
                                      n_target_genes = n_target_genes,
                                      n_regulator_genes = n_regulator_genes,
                                      n_cells = n_cells,
                                      n_cells_test = n_cells_test,
                                      regulator_means = regulator_means,
                                      coefficient_means = coefficient_means,
                                      coefficient_sds = coefficient_sds,
                                      disturbed_fraction = disturbed_fraction,
                                      testing_penalization_data_gen = testing_penalization_data_gen,
                                      generated_data = generated_data,
                                      ind_targetgenes = ind_targetgenes,
                                      ind_reggenes = ind_reggenes,
                                      disturbed_initial_cell_clust = disturbed_initial_cell_clust,
                                      disturbed_initial_cell_clust_test = disturbed_initial_cell_clust_test,
                                      biclust_input_data = biclust_input_data,
                                      biclust_input_data_test = biclust_input_data_test,
                                      n_total_cells = n_total_cells,
                                      n_total_cells_test = n_total_cells_test,
                                      cell_id = cell_id,
                                      cell_id_test = cell_id_test,
                                      correct_clustering = correct_clustering,
                                      correct_clustering_test = correct_clustering_test,
                                      trivial_regulator_networks = trivial_regulator_networks,
                                      seed = seed)))
    return(scenario_list)
  }



  generated_data <- generate_dummy_data_for_cell_clustering(
    n_cell_clusters = n_cell_clusters,
    n_target_gene_clusters = n_target_gene_clusters,
    # Number of target gene clusters in each cell cluster
    n_target_genes = n_target_genes,
    #from vignette
    n_regulator_genes = n_regulator_genes,
    # from vignette
    n_cells = n_cells,
    regulator_means = regulator_means,
    # For generating dummy data, regulator mean in each cell cluster
    coefficient_means <- coefficient_means,
    # For generating dummy data, coefficient means in each cell cluster
    coefficient_sds <- coefficient_sds,
    disturbed_fraction = disturbed_fraction,
    # Value between 0 and 1. How large portion of cells should move to other cell clusters.
    plot_stuff = FALSE,
    plot_suffix = plot_suffix,
    testing_penalization = testing_penalization_data_gen,
    generate_counts             = FALSE,
    check_results               = FALSE,
    trivial_regulator_networks  = trivial_regulator_networks,
    pearson_regulators          = TRUE
  )


  # Set some more variables -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # Because "dat <- cbind(Z_t, Z_r)" in generate_dummy_data_for_cell_clustering
  ind_targetgenes <- which(c(rep(1, n_target_genes), rep(0, n_regulator_genes)) == 1)
  ind_reggenes <- which(c(rep(0, n_target_genes), rep(1, n_regulator_genes)) == 1)

  disturbed_initial_cell_clust <- factor(generated_data$disturbed_initial_cell_clust)

  biclust_input_data <- generated_data$dat
  colnames(biclust_input_data) <- c(paste0("t", 1:n_target_genes),
                                    paste0("r", 1:n_regulator_genes))
  biclust_input_data <- tibble::as_tibble(biclust_input_data)

  n_total_cells <- sum(n_cells)
  cell_id <- 1:n_total_cells


  # Construct correct biclust matrix --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  true_biclust_matrix <- construct_biclust_matrix(cell_cluster_allocation = generated_data$true_cell_clust,
                                                  target_gene_allocation= generated_data$true_target_gene_allocation)
  correct_clustering <- as.vector(true_biclust_matrix)


  # Plot heatmap of data --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  if(plot_heatmap){
    n <- length(unique(as.vector(true_biclust_matrix)))
    if(n!=max(as.vector(true_biclust_matrix))){
      print("Warning")
    }
    regions <- seq(1, n, length.out = n + 1)
    middle_of_regions <- (regions[-1] + regions[-length(regions)]) / 2
    odd_number_larger <- ifelse(n %% 2 == 0, n + 1, n)
    if(n %% 2 == 1){
      keep_these_colors = 1:n
    }else{
      keep_these_colors <- setdiff(1:odd_number_larger, (odd_number_larger + 1) / 2 + 1)
    }
    constructed_plot <- rasterVis::levelplot(true_biclust_matrix, att = n,
                                             # col.regions = rainbow(odd_number_larger),
                                             colorkey = list(at = regions,
                                                             # col=rainbow(odd_number_larger)[keep_these_colors],
                                                             labels = list(at = middle_of_regions, labels = as.character(1:n))),
                                             xlab = 'Cells',
                                             ylab = 'Target genes',
                                             main='Generated data')

    pdf(file.path(output_path, paste0("heatmap_generateddata_",plot_suffix,".pdf")), width = 12, height = 6)
    print(constructed_plot)
    dev.off()

  }

  st <- split_traintest( generated_data,
                         disturbed_initial_cell_clust,
                         biclust_input_data,
                         n_cells,
                         n_total_cells,  # num 400
                         cell_id,  # num [1:400] 1 2 3 4 5 6 7 8 9 10 ...
                         correct_clustering)
  generated_data <- st$generated_data
  disturbed_initial_cell_clust <- st$disturbed_initial_cell_clust
  disturbed_initial_cell_clust_test <- st$disturbed_initial_cell_clust_test
  biclust_input_data <- st$biclust_input_data
  biclust_input_data_test <- st$biclust_input_data_test
  n_total_cells <- st$n_total_cells
  n_total_cells_test <- st$n_total_cells_test
  cell_id <- st$cell_id
  cell_id_test <- st$cell_id_test
  correct_clustering <- st$correct_clustering
  correct_clustering_test <- st$correct_clustering_test
  n_cells <- st$n_cells
  n_cells_test <- st$n_cells_test

  scenario_list <- add_scenario_data(scenario_list,
                                     description, # : chr "Simple"
                                     n_cell_clusters,# : int 2
                                     n_target_gene_clusters, # : int [1:2] 4 2
                                     n_target_genes, # : int 100
                                     n_regulator_genes, #  : int 6
                                     n_cells, #  : num [1:2] 200 200
                                     n_cells_test, #  : num [1:2] 200 200
                                     regulator_means, #  : num [1:2] 0 0
                                     coefficient_means, #  :List of 2
                                     # ..$ : num [1:4] 1 3 5 7
                                     # ..$ : num [1:2] 10 20
                                     coefficient_sds, #  :List of 2
                                     # ..$ : num [1:4] 0.01 0.01 0.01 0.01
                                     # ..$ : num [1:2] 0.01 0.01
                                     disturbed_fraction,  # num [1:2] 0.1 0.5
                                     testing_penalization_data_gen, # num [1:2] 0.1 0.5
                                     generated_data, #  generated_data               :List of 9
                                     # ..$ disturbed_initial_cell_clust: int [1:400] 1 1 1 1 1 1 1 1 1 1 ...
                                     # ..$ initial_cell_clust          : int [1:400] 1 1 1 1 1 1 1 1 1 1 ...
                                     # ..$ true_cell_clust             : int [1:400] 1 1 1 1 1 1 1 1 1 1 ...
                                     # ..$ true_betas                  :List of 2
                                     # .. ..$ :List of 4
                                     # .. .. ..$ : num [1:6, 1:100] 0.994 0 0 0 0 ...
                                     # .. .. ..$ : num [1:6, 1:100] 0 0 0 0 0 0 0 0 0 0 ...
                                     # .. .. ..$ : num [1:6, 1:100] 0 0 0 0 0 0 0 0 0 0 ...
                                     # .. .. ..$ : num [1:6, 1:100] 0 0 0 0 0 0 0 0 0 0 ...
                                     # .. ..$ :List of 2
                                     # .. .. ..$ : num [1:6, 1:100] 0 0 0 0 9.99 ...
                                     # .. .. ..$ : num [1:6, 1:100] 0 0 0 0 0 0 0 0 0 0 ...
                                     # ..$ dat                         : num [1:400, 1:106] -2.19 -3.14 1.58 -2.82 -2.51 ...
                                     # ..$ true_target_gene_allocation :List of 2
                                     # .. ..$ : int [1:100] 1 1 1 1 1 1 1 1 1 1 ...
                                     # .. ..$ : int [1:100] 1 1 1 1 1 1 1 1 1 1 ...
                                     # ..$ true_Pi                     :List of 2
                                     # .. ..$ : num [1:4, 1:100] 1 0 0 0 1 0 0 0 1 0 ...
                                     # .. ..$ : num [1:2, 1:100] 1 0 1 0 1 0 1 0 1 0 ...
                                     # ..$ true_S                      :List of 2
                                     # .. ..$ : num [1:4, 1:6] 1 0 0 0 0 1 0 0 0 0 ...
                                     # .. ..$ : num [1:2, 1:6] 0 0 0 0 0 0 0 0 1 0 ...
                                     # ..$ counts                      : logi NA
                                     ind_targetgenes, # int [1:100] 1 2 3 4 5 6 7 8 9 10 ...
                                     ind_reggenes, # int [1:6] 101 102 103 104 105 106
                                     disturbed_initial_cell_clust, # Factor w/ 2 levels "1","2": 1 1 1 1 1 1 1 1 1 1 ...
                                     disturbed_initial_cell_clust_test, # Factor w/ 2 levels "1","2": 1 1 1 1 1 1 1 1 1 1 ...
                                     biclust_input_data,  # tibble [400 × 106] (S3: tbl_df/tbl/data.frame)
                                     biclust_input_data_test,  # tibble [400 × 106] (S3: tbl_df/tbl/data.frame)
                                     n_total_cells,  # num 400
                                     n_total_cells_test,  # num 400
                                     cell_id,  # num [1:400] 1 2 3 4 5 6 7 8 9 10 ...
                                     cell_id_test,  # num [1:400] 401 402 403 404 405 406 407 408 409 410 ...
                                     correct_clustering,  # biclust clustering: int [1:40000] 1 1 1 1 1 1 1 1 1 1 ...
                                     correct_clustering_test,  # biclust clustering: int [1:40000] 1 1 1 1 1 1 1 1 1 1 ...
                                     trivial_regulator_networks, # logi TRUE
                                     seed) # num 1234



  return(scenario_list)
}



# Create an empty list to store all scenarios ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if (!file.exists(file.path(output_path_rds, "sim_data.rds")) | redo_flag) {



  scenarios <- list()

  # Generate one simple scenario ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  scenarios <- create_scenario(scenario_list = scenarios,
                               description = "Simple",
                               n_cell_clusters = 2,
                               n_target_gene_clusters = c(4, 2),  # Number of target gene clusters in each cell cluster
                               n_target_genes = 100,
                               n_regulator_genes = 6,
                               n_cells = c(200, 200)*2,
                               regulator_means = c(0, 0), # For generating dummy data, regulator mean in each cell cluster
                               coefficient_means = list(c(1, 3, 5, 7), c(10, 20)),  # For generating dummy data, coefficient means in each cell cluster
                               coefficient_sds = list(c(0.01, 0.01, 0.01, 0.01), c(0.01, 0.01)),
                               disturbed_fraction = 0.1,  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
                               testing_penalization_data_gen = c(0.1, 0.5),
                               trivial_regulator_networks = TRUE,
                               seed = 1234,
                               plot_heatmap=TRUE,
                               plot_suffix="simple-scenario-1")

  scenarios <- create_scenario(scenario_list = scenarios,
                               description = "Simple",
                               n_cell_clusters = 3,
                               n_target_gene_clusters = c(4, 2, 2),  # Number of target gene clusters in each cell cluster
                               n_target_genes = 120,
                               n_regulator_genes = 8,
                               n_cells = c(220, 220, 220)*2,
                               regulator_means = c(0, 0, 0), # For generating dummy data, regulator mean in each cell cluster
                               coefficient_means = list(c(1, 3, 5, 7), c(10, 12), c(22, 24)),  # For generating dummy data, coefficient means in each cell cluster
                               coefficient_sds = list(c(0.01, 0.01, 0.01, 0.01), c(0.01, 0.01), c(0.01, 0.01)),
                               disturbed_fraction = 0.1,  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
                               testing_penalization_data_gen = c(0.5, 0.5, 0.5),
                               trivial_regulator_networks = TRUE,
                               seed = 1234)

  scenarios <- create_scenario(scenario_list = scenarios,
                               description = "Simple",
                               n_cell_clusters = 4,
                               n_target_gene_clusters = c(4, 2, 2, 3),  # Number of target gene clusters in each cell cluster
                               n_target_genes = 120,
                               n_regulator_genes = 11,
                               n_cells = c(220, 220, 220, 220)*2,
                               regulator_means = c(0, 0, 0, 0), # For generating dummy data, regulator mean in each cell cluster
                               coefficient_means = list(c(1, 3, 5, 7), c(10, 12), c(22, 24), c(33, 35, 37)),  # For generating dummy data, coefficient means in each cell cluster
                               coefficient_sds = list(c(0.01, 0.01, 0.01, 0.01), c(0.01, 0.01), c(0.01, 0.01), c(0.01, 0.01, 0.01)),
                               disturbed_fraction = 0.1,  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
                               testing_penalization_data_gen = c(0.5, 0.5, 0.5, 0.5),
                               trivial_regulator_networks = TRUE,
                               seed = 1234)

  scenarios <- create_scenario(scenario_list = scenarios,
                               description = "Simple",
                               n_cell_clusters = 4,
                               n_target_gene_clusters = c(4, 3, 2, 3),  # Number of target gene clusters in each cell cluster
                               n_target_genes = 160,
                               n_regulator_genes = 12,
                               n_cells = c(220, 220, 220, 220)*2,
                               regulator_means = c(0, 0, 0, 0), # For generating dummy data, regulator mean in each cell cluster
                               coefficient_means = list(c(1, 3, 5, 7), c(15, 17, 13), c(22, 24), c(33, 35, 37)),  # For generating dummy data, coefficient means in each cell cluster
                               coefficient_sds = list(c(0.01, 0.01, 0.01, 0.01), c(0.01, 0.01, 0.01), c(0.01, 0.01), c(0.01, 0.01, 0.01)),
                               disturbed_fraction = 0.1,  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
                               testing_penalization_data_gen = c(0.5, 0.5, 0.5, 0.5),
                               trivial_regulator_networks = TRUE,
                               seed = 1234)

  scenarios <- create_scenario(scenario_list = scenarios,
                               description = "Simple",
                               n_cell_clusters = 4,
                               n_target_gene_clusters = c(2, 2, 2, 5),  # Number of target gene clusters in each cell cluster
                               n_target_genes = 240,
                               n_regulator_genes = 11,
                               n_cells = c(250, 250, 250, 250)*2,
                               regulator_means = c(0, 0, 0, 0), # For generating dummy data, regulator mean in each cell cluster
                               coefficient_means = list(c(1, 3), c(15, 17), c(22, 24), c(33, 35, 37, 39, 41)),  # For generating dummy data, coefficient means in each cell cluster
                               coefficient_sds = list(c(0.01, 0.01), c(0.01, 0.01), c(0.01, 0.01), c(0.01, 0.01, 0.01, 0.01, 0.01)),
                               disturbed_fraction = 0.1,  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
                               testing_penalization_data_gen = c(0.5, 0.5, 0.5, 0.5),
                               trivial_regulator_networks = TRUE,
                               seed = 1234)

  scenarios <- create_scenario(scenario_list = scenarios,
                               description = "Simple",
                               n_cell_clusters = 2,
                               n_target_gene_clusters = c(4, 2),  # Number of target gene clusters in each cell cluster
                               n_target_genes = 100,
                               n_regulator_genes = 6,
                               n_cells = c(200, 200)*2,
                               regulator_means = c(0, 0), # For generating dummy data, regulator mean in each cell cluster
                               coefficient_means = list(c(1, 3, 5, 7), c(10, 20)),  # For generating dummy data, coefficient means in each cell cluster
                               coefficient_sds = list(c(0.06, 0.03, 0.04, 0.02), c(0.05, 0.02)),
                               disturbed_fraction = 0.1,  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
                               testing_penalization_data_gen = c(0.1, 0.5),
                               trivial_regulator_networks = TRUE,
                               seed = 1234)

  scenarios <- create_scenario(scenario_list = scenarios,
                               description = "Simple",
                               n_cell_clusters = 2,
                               n_target_gene_clusters = c(4, 2),  # Number of target gene clusters in each cell cluster
                               n_target_genes = 100,
                               n_regulator_genes = 6,
                               n_cells = c(400, 200)*2,
                               regulator_means = c(0, 0), # For generating dummy data, regulator mean in each cell cluster
                               coefficient_means = list(c(1, 3, 5, 7), c(10, 20)),  # For generating dummy data, coefficient means in each cell cluster
                               coefficient_sds = list(c(0.1, 0.03, 0.04, 0.02), c(0.05, 0.02)),
                               disturbed_fraction = 0.1,  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
                               testing_penalization_data_gen = c(0.1, 0.5),
                               trivial_regulator_networks = TRUE,
                               seed = 1234)

  scenarios <- create_scenario(scenario_list = scenarios,
                               description = "Simple",
                               n_cell_clusters = 2,
                               n_target_gene_clusters = c(7, 3),  # Number of target gene clusters in each cell cluster
                               n_target_genes = 100,
                               n_regulator_genes = 10,
                               n_cells = c(600, 300)*2,
                               regulator_means = c(0, 0), # For generating dummy data, regulator mean in each cell cluster
                               coefficient_means = list(c(1, 3, 5, 7, 9, 11, 13), c(22, 20, 24)),  # For generating dummy data, coefficient means in each cell cluster
                               coefficient_sds = list(c(0.1, 0.03, 0.04, 0.02, 0.01, 0.01, 0.01), c(0.05, 0.02, 0.01)),
                               disturbed_fraction = 0.1,  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
                               testing_penalization_data_gen = c(0.1, 0.5),
                               trivial_regulator_networks = TRUE,
                               seed = 1234)

  scenarios <- create_scenario(scenario_list = scenarios,
                               description = "Simple",
                               n_cell_clusters = 3,
                               n_target_gene_clusters = c(2, 2, 2),  # Number of target gene clusters in each cell cluster
                               n_target_genes = 120,
                               n_regulator_genes = 6,
                               n_cells = c(300, 220, 260)*2,
                               regulator_means = c(0, 0, 0), # For generating dummy data, regulator mean in each cell cluster
                               coefficient_means = list(c(1, 3), c(10, 12), c(22, 24)),  # For generating dummy data, coefficient means in each cell cluster
                               coefficient_sds = list(c(0.02, 0.06), c(0.04, 0.03), c(0.02, 0.04)),
                               disturbed_fraction = 0.1,  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
                               testing_penalization_data_gen = c(0.5, 0.5, 0.5),
                               trivial_regulator_networks = TRUE,
                               seed = 1234)

  scenarios <- create_scenario(scenario_list = scenarios,
                               description = "Simple",
                               n_cell_clusters = 3,
                               n_target_gene_clusters = c(2, 10, 2),  # Number of target gene clusters in each cell cluster
                               n_target_genes = 220,
                               n_regulator_genes = 14,
                               n_cells = c(300, 2200, 260)*2,
                               regulator_means = c(0, 0, 0), # For generating dummy data, regulator mean in each cell cluster
                               coefficient_means = list(c(1, 3), c(10, 12, 14, 16, 18, 20, 22, 24, 26, 28), c(51, 53)),  # For generating dummy data, coefficient means in each cell cluster
                               coefficient_sds = list(c(0.02, 0.06), c(0.04, 0.03, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01), c(0.02, 0.04)),
                               disturbed_fraction = 0.1,  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
                               testing_penalization_data_gen = c(0.5, 0.5, 0.5),
                               trivial_regulator_networks = TRUE,
                               seed = 1234)



  # Generate multiple complicated examples --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  scenarios <- create_scenario(scenario_list = scenarios,
                               description = "Complex",
                               n_cell_clusters = 2,
                               n_target_gene_clusters = c(3, 2),  # Number of target gene clusters in each cell cluster
                               n_target_genes = 300,
                               n_regulator_genes = 20,
                               n_cells = c(600, 400)*2,
                               regulator_means = c(0, 0), # For generating dummy data, regulator mean in each cell cluster
                               coefficient_means = list(c(1, 2, 3), c(18, 20)),  # For generating dummy data, coefficient means in each cell cluster
                               coefficient_sds = list(c(0.01, 0.01, 0.01), c(0.03, 0.01)),
                               disturbed_fraction = 0.1,  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
                               testing_penalization_data_gen = c(0.5, 0.5),
                               trivial_regulator_networks = FALSE,
                               seed = 1234)

  scenarios <- create_scenario(scenario_list = scenarios,
                               description = "Complex",
                               n_cell_clusters = 2,
                               n_target_gene_clusters = c(4, 2),  # Number of target gene clusters in each cell cluster
                               n_target_genes = 300,
                               n_regulator_genes = 20,
                               n_cells = c(500, 300)*2,
                               regulator_means = c(0, 0), # For generating dummy data, regulator mean in each cell cluster
                               coefficient_means = list(c(1, 3, 5, 7), c(18, 20)),  # For generating dummy data, coefficient means in each cell cluster
                               coefficient_sds = list(c(0.01, 0.01, 0.01, 0.01), c(0.01, 0.01)),
                               disturbed_fraction = 0.1,  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
                               testing_penalization_data_gen = c(0.5, 0.5),
                               trivial_regulator_networks = FALSE,
                               seed = 123)

  scenarios <- create_scenario(scenario_list = scenarios,
                               description = "Complex",
                               n_cell_clusters = 2,
                               n_target_gene_clusters = c(3, 2),  # Number of target gene clusters in each cell cluster
                               n_target_genes = 300,
                               n_regulator_genes = 20,
                               n_cells = c(600, 400)*2,
                               regulator_means = c(0, 0), # For generating dummy data, regulator mean in each cell cluster
                               coefficient_means = list(c(1, 2, 3), c(18, 20)),  # For generating dummy data, coefficient means in each cell cluster
                               coefficient_sds = list(c(0.01, 0.01, 0.01), c(0.03, 0.01)),
                               disturbed_fraction = 0.3,  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
                               testing_penalization_data_gen = c(0.5, 0.5),
                               trivial_regulator_networks = FALSE,
                               seed = 124)

  scenarios <- create_scenario(scenario_list = scenarios,
                               description = "Complex",
                               n_cell_clusters = 3,
                               n_target_gene_clusters = c(3, 3, 3),  # Number of target gene clusters in each cell cluster
                               n_target_genes = 300,
                               n_regulator_genes = 20,
                               n_cells = c(500, 300, 500)*2,
                               regulator_means = c(0, 0, 0), # For generating dummy data, regulator mean in each cell cluster
                               coefficient_means = list(c(1, 2, 3), c(22, 20, 24), c(10,12,14)),  # For generating dummy data, coefficient means in each cell cluster
                               coefficient_sds = list(c(0.03, 0.01, 0.03), c(0.02, 0.04, 0.01), c(0.02, 0.04, 0.01)),
                               disturbed_fraction = 0.5,  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
                               testing_penalization_data_gen = c(0.5, 0.5, 0.5),
                               trivial_regulator_networks = FALSE,
                               seed = 1236)

  scenarios <- create_scenario(scenario_list = scenarios,
                               description = "Complex",
                               n_cell_clusters = 4,
                               n_target_gene_clusters = c(3, 3, 3, 5),  # Number of target gene clusters in each cell cluster
                               n_target_genes = 500,
                               n_regulator_genes = 30,
                               n_cells = c(500, 300, 500, 300)*2,
                               regulator_means = c(0, 0, 0, 0), # For generating dummy data, regulator mean in each cell cluster
                               coefficient_means = list(c(1, 2, 3), c(22, 20, 24), c(10,12,14), c(40,42,44,46,48)),  # For generating dummy data, coefficient means in each cell cluster
                               coefficient_sds = list(c(0.03, 0.01, 0.03), c(0.02, 0.04, 0.01), c(0.02, 0.04, 0.01), c(0.01, 0.02, 0.012, 0.011, 0.005)),
                               disturbed_fraction = 0.101,  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
                               testing_penalization_data_gen = c(0.5, 0.5, 0.5, 0.5),
                               trivial_regulator_networks = FALSE,
                               seed = 126)

  scenarios <- create_scenario(scenario_list = scenarios,
                               description = "Complex",
                               n_cell_clusters = 2,
                               n_target_gene_clusters = c(2, 2),  # Number of target gene clusters in each cell cluster
                               n_target_genes = 100,
                               n_regulator_genes = 8,
                               n_cells = c(340, 340)*2,
                               regulator_means = c(0, 0), # For generating dummy data, regulator mean in each cell cluster
                               coefficient_means = list(c(5, 3), c(18, 20)),  # For generating dummy data, coefficient means in each cell cluster
                               coefficient_sds = list(c(0.01, 0.01), c(0.03, 0.05)),
                               disturbed_fraction = 0.23,  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
                               testing_penalization_data_gen = c(0.5, 0.5),
                               trivial_regulator_networks = FALSE,
                               seed = 1224)

  scenarios <- create_scenario(scenario_list = scenarios,
                               description = "Complex",
                               n_cell_clusters = 2,
                               n_target_gene_clusters = c(4, 4),  # Number of target gene clusters in each cell cluster
                               n_target_genes = 120,
                               n_regulator_genes = 12,
                               n_cells = c(440, 440)*2,
                               regulator_means = c(0, 0), # For generating dummy data, regulator mean in each cell cluster
                               coefficient_means = list(c(15, 13, 14, 10), c(31, 30, 28, 29)),  # For generating dummy data, coefficient means in each cell cluster
                               coefficient_sds = list(c(0.03, 0.012, 0.011, 0.015), c(0.03, 0.05, 0.01, 0.01)),
                               disturbed_fraction = 0.33,  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
                               testing_penalization_data_gen = c(0.5, 0.5),
                               trivial_regulator_networks = FALSE,
                               seed = 124)

  scenarios <- create_scenario(scenario_list = scenarios,
                               description = "Complex",
                               n_cell_clusters = 2,
                               n_target_gene_clusters = c(2, 3),  # Number of target gene clusters in each cell cluster
                               n_target_genes = 130,
                               n_regulator_genes = 15,
                               n_cells = c(500, 500)*2,
                               regulator_means = c(0, 0), # For generating dummy data, regulator mean in each cell cluster
                               coefficient_means = list(c(15, 13), c(30, 28, 29)),  # For generating dummy data, coefficient means in each cell cluster
                               coefficient_sds = list(c(0.03, 0.012), c(0.03, 0.05, 0.01)),
                               disturbed_fraction = 0.23,  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
                               testing_penalization_data_gen = c(0.5, 0.5),
                               trivial_regulator_networks = FALSE,
                               seed = 1214)

  scenarios <- create_scenario(scenario_list = scenarios,
                               description = "Complex",
                               n_cell_clusters = 4,
                               n_target_gene_clusters = c(2, 3, 2, 2),  # Number of target gene clusters in each cell cluster
                               n_target_genes = 130,
                               n_regulator_genes = 15,
                               n_cells = c(500, 500, 400, 400)*2,
                               regulator_means = c(0, 0, 0, 0), # For generating dummy data, regulator mean in each cell cluster
                               coefficient_means = list(c(15, 13), c(30, 28, 29), c(43, 45), c(20,22)),  # For generating dummy data, coefficient means in each cell cluster
                               coefficient_sds = list(c(0.03, 0.012), c(0.03, 0.05, 0.01), c(0.01, 0.01), c(0.01, 0.01)),
                               disturbed_fraction = 0.13,  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
                               testing_penalization_data_gen = c(0.5, 0.5, 0.5, 0.5),
                               trivial_regulator_networks = FALSE,
                               seed = 1241)

  scenarios <- create_scenario(scenario_list = scenarios,
                               description = "Complex",
                               n_cell_clusters = 2,
                               n_target_gene_clusters = c(4, 2),  # Number of target gene clusters in each cell cluster
                               n_target_genes = 100,
                               n_regulator_genes = 8,
                               n_cells = c(220, 220)*2,
                               regulator_means = c(0, 0), # For generating dummy data, regulator mean in each cell cluster
                               coefficient_means = list(c(1, 3, 5, 7), c(25, 27)),  # For generating dummy data, coefficient means in each cell cluster
                               coefficient_sds = list(c(0.03, 0.02, 0.01, 0.01), c(0.01, 0.01)),
                               disturbed_fraction = 0.56,  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
                               testing_penalization_data_gen = c(0.5, 0.5),
                               trivial_regulator_networks = FALSE,
                               seed = 631)

  saveRDS(scenarios, file.path(output_path_rds, "sim_data.rds"))

} else {


  scenarios <- readRDS(file.path(output_path_rds, "sim_data.rds"))

}






#
# for(iter in 1:20){
#   n_cell_clusters               <- sample(2:10, 1)
#   n_target_genes                <- 100
#
#   # Number of target gene clusters in each cell cluster
#   n_target_gene_clusters        <- sample(2:(n_target_genes/10), n_cell_clusters, replace = TRUE)
#
#   coefficient_sds               <- lapply(n_target_gene_clusters, FUN = function(cellClust){
#     # sample(seq(from= 0.01, to= 1, length.out = 99),
#     #        cellClust, replace = T)
#
#     # instead maybe sample a value for each cell cluster,
#     # and draw some different means around those values
#     rep(sample(seq(from= 0.01, to= 1, length.out = 99), 1),
#         cellClust)
#
#   }
#   )
#
#   # For generating dummy data, coefficient means in each cell cluster
#   coefficient_means <- lapply(n_target_gene_clusters, FUN = function(cellClust){
#     #         sample(1:40,cellClust, replace = T)
#     # instead maybe sample a value for each cell cluster,
#     # and draw some different means around those values
#     rep(sample(seq(from= -10, to= 10, length.out = 99), 1),
#         cellClust)
#     #maybe also perturb each mean so that hey are separated by at least the corresponding sd
#   }
#   )
#
#
#   scenarios <- create_scenario(scenario_list = scenarios,
#                                description = "Complex",
#                                n_cell_clusters = n_cell_clusters,
#                                n_target_gene_clusters = n_target_gene_clusters,
#                                n_target_genes = n_target_genes,
#                                n_regulator_genes = 10,
#                                n_cells = sample(100:1000, n_cell_clusters, replace = TRUE),
#                                regulator_means = rep(0, n_cell_clusters), # For generating dummy data, regulator mean in each cell cluster
#                                coefficient_means = coefficient_means,
#                                coefficient_sds = coefficient_sds,
#                                disturbed_fraction = 0.1,  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
#                                testing_penalization_data_gen = rep(0.5, n_cell_clusters),  # one per cell cluster
#                                trivial_regulator_networks = FALSE,
#                                seed = 1234)
#
# }
#




