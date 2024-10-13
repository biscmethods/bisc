# !/usr/bin/Rscript
rm(list = ls())

library(patchwork)
library(rasterVis)
library(cowplot)
sink()

# # options(warn=2)  # To convert warning messages into error messages which display row of error. For debugging.
#
# # Get absolute path where script is located, by using relative paths.
# demo_path <- here::here("demo")
R_path <- here::here("R")
output_path <- here::here("demo/simulated_data_to_stats_pipeline")
# path_data <- here::here('data')
#


source(file.path(R_path, "generate_dummy_data_for_cell_clustering.R"))



# Function that creates data for a scenario and returns a list of all variables needed for both generating the same data and run analysis -------------------------------------------------------------------------------------------------
create_scenario <- function(scenario_list,
                            description,
                            plot_heatmap = FALSE,
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
                                regulator_means,
                                coefficient_means,
                                coefficient_sds,
                                disturbed_fraction,
                                testing_penalization_data_gen,
                                generated_data,
                                ind_targetgenes,
                                ind_reggenes,
                                disturbed_initial_cell_clust,
                                biclust_input_data,
                                n_total_cells,
                                cell_id,
                                correct_clustering,
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
                                      regulator_means = regulator_means,
                                      coefficient_means = coefficient_means,
                                      coefficient_sds = coefficient_sds,
                                      disturbed_fraction = disturbed_fraction,
                                      testing_penalization_data_gen = testing_penalization_data_gen,
                                      generated_data = generated_data,
                                      ind_targetgenes = ind_targetgenes,
                                      ind_reggenes = ind_reggenes,
                                      disturbed_initial_cell_clust = disturbed_initial_cell_clust,
                                      biclust_input_data = biclust_input_data,
                                      n_total_cells = n_total_cells,
                                      cell_id = cell_id,
                                      correct_clustering = correct_clustering,
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
    plot_suffix = "",
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
  cell_cluster_allocation <- generated_data$true_cell_clust
  target_gene_allocation <- generated_data$true_target_gene_allocation
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
    print(constructed_plot)
  }

  scenario_list <- add_scenario_data(scenario_list,
                                     description,
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
                                     generated_data,
                                     ind_targetgenes,
                                     ind_reggenes,
                                     disturbed_initial_cell_clust,
                                     biclust_input_data,
                                     n_total_cells,
                                     cell_id,
                                     correct_clustering,
                                     trivial_regulator_networks,
                                     seed)
  return(scenario_list)
}



# Create an empty list to store all scenarios ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

scenarios <- list()

# Generate one simple scenario ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

scenarios <- create_scenario(scenario_list = scenarios,
                             description = "Simple",
                             n_cell_clusters = 2,
                             n_target_gene_clusters = c(4, 2),  # Number of target gene clusters in each cell cluster
                             n_target_genes = 100,
                             n_regulator_genes = 6,
                             n_cells = c(200, 200),
                             regulator_means = c(0, 0), # For generating dummy data, regulator mean in each cell cluster
                             coefficient_means = list(c(1, 3, 5, 7), c(10, 20)),  # For generating dummy data, coefficient means in each cell cluster
                             coefficient_sds = list(c(0.01, 0.01, 0.01, 0.01), c(0.01, 0.01)),
                             disturbed_fraction = 0.1,  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
                             testing_penalization_data_gen = c(0.1, 0.5),
                             trivial_regulator_networks = TRUE,
                             seed = 1234)

scenarios <- create_scenario(scenario_list = scenarios,
                             description = "Simple",
                             n_cell_clusters = 3,
                             n_target_gene_clusters = c(4, 2, 2),  # Number of target gene clusters in each cell cluster
                             n_target_genes = 120,
                             n_regulator_genes = 8,
                             n_cells = c(220, 220, 220),
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
                             n_cells = c(220, 220, 220, 220),
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
                             n_cells = c(220, 220, 220, 220),
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
                             n_cells = c(250, 250, 250, 250),
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
                             n_cells = c(200, 200),
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
                             n_cells = c(400, 200),
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
                             n_cells = c(600, 300),
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
                             n_cells = c(300, 220, 260),
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
                             n_cells = c(300, 2200, 260),
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
                             n_cells = c(600, 400),
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
                             n_cells = c(500, 300),
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
                             n_cells = c(600, 400),
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
                             n_cells = c(500, 300, 500),
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
                             n_cells = c(500, 300, 500, 300),
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
                             n_cells = c(340, 340),
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
                             n_cells = c(440, 440),
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
                             n_cells = c(500, 500),
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
                             n_cells = c(500, 500, 400, 400),
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
                             n_cells = c(220, 220),
                             regulator_means = c(0, 0), # For generating dummy data, regulator mean in each cell cluster
                             coefficient_means = list(c(1, 3, 5, 7), c(25, 27)),  # For generating dummy data, coefficient means in each cell cluster
                             coefficient_sds = list(c(0.03, 0.02, 0.01, 0.01), c(0.01, 0.01)),
                             disturbed_fraction = 0.56,  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
                             testing_penalization_data_gen = c(0.5, 0.5),
                             trivial_regulator_networks = FALSE,
                             seed = 631)
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

