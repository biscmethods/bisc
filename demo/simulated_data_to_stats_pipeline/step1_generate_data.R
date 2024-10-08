# #!/usr/bin/Rscript
# rm(list = ls())
#

library(patchwork)
library(rasterVis)
library(cowplot)
# sink()
#
# # options(warn=2)  # To convert warning messages into error messages which display row of error. For debugging.
#
# # Get absolute path where script is located, by using relative paths.
# demo_path <- here::here("demo")
# R_path <- here::here("R")
# output_path <- here::here("demo/simulated_data_to_stats_pipeline")
# path_data <- here::here('data')
#
#
# redo_flag <- TRUE

source(file.path(R_path, "generate_dummy_data_for_cell_clustering.R"))


################################################################
############ Generate initial simple example ###################
################################################################



# Set seed for example
set.seed(1234)

# Set variables ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
n_cell_clusters <- 2
n_target_gene_clusters <- c(4, 2)  # Number of target gene clusters in each cell cluster
n_target_genes <- 100
n_regulator_genes <- 6
n_cells <- c(200, 200)
regulator_means <- c(0, 0) # For generating dummy data, regulator mean in each cell cluster
coefficient_means <- list(c(1, 3, 5, 7), c(10, 20))  # For generating dummy data, coefficient means in each cell cluster
coefficient_sds <- list(c(0.01, 0.01, 0.01, 0.01), c(0.01, 0.01))
disturbed_fraction <- 0.1  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
testing_penalization_data_gen <- c(0.1, 0.5)

if (!file.exists(file.path(output_path, "env_sim_simple_data_biclust_sc.rds")) |
    redo_flag) {
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
    plot_suffix = "Simple",
    testing_penalization = testing_penalization_data_gen,
    generate_counts             = FALSE,
    check_results               = FALSE,
    trivial_regulator_networks  = TRUE,
    pearson_regulators          = TRUE
  )


  saveRDS(generated_data,
          file.path(output_path, "env_sim_simple_data_biclust_sc.rds"))

} else {
  generated_data <- readRDS(file.path(output_path, "env_sim_simple_data_biclust_sc.rds"))

}


# Plot Counts -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# devtools::install_github("jokergoo/ComplexHeatmap")
# d <- generated_data$counts
# d_u <- sort(unique(as.vector(d)))
# colors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlBu"))(max(length(d_u), max(d_u)))
# colors[1] <- "#000000"
# colors <- structure(colors, names = as.character(d_u))
# ComplexHeatmap::Heatmap(d,
#                         name = "Count",
#                         col = colors,
#                         column_title = "Generated data counts",
#                         cluster_rows = FALSE,
#                         cluster_columns = FALSE)
#--------------

# Because "dat <- cbind(Z_t, Z_r)" in generate_dummy_data_for_cell_clustering
ind_targetgenes <- which(c(rep(1, n_target_genes), rep(0, n_regulator_genes)) == 1)
ind_reggenes <- which(c(rep(0, n_target_genes), rep(1, n_regulator_genes)) == 1)


disturbed_initial_cell_clust <- factor(generated_data$disturbed_initial_cell_clust)

biclust_input_data <- generated_data$dat
colnames(biclust_input_data) <- c(paste0("t", 1:n_target_genes),
                                  paste0("r", 1:n_regulator_genes))
biclust_input_data <- tibble::as_tibble(biclust_input_data)

# Set up some variables
n_cell_clusters <- length(unique(disturbed_initial_cell_clust))
n_target_genes <- length(ind_targetgenes)
n_regulator_genes <- length(ind_reggenes)
n_total_cells <- sum(n_cells)
cell_id <- 1:n_total_cells


# Construct heatmap for generated data ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

cell_cluster_allocation <- generated_data$true_cell_clust
target_gene_allocation <- generated_data$true_target_gene_allocation

result_matrix <- matrix(0, nrow = n_total_cells, ncol = n_target_genes)

# Check each entry in the matrix (each cell and gene pair), and assign a string-number to each unique "cell-cluster-gene-module".
for (i in 1:n_total_cells) {
  cluster <- cell_cluster_allocation[i]
  gene_allocation <- target_gene_allocation[[cluster]][1:n_target_genes]
  gene_allocation[gene_allocation==-1] <- 0
  # Create unique numbers for each pair
  result_matrix[i, ] <- paste0(cluster, gene_allocation)
}

# Convert the string-numbers to numeric matrix (starting from 1 this time)
result_matrix <- matrix(as.numeric(result_matrix), nrow = n_total_cells, ncol = n_target_genes)
result_matrix <- matrix(as.integer(as.factor(result_matrix)), nrow=nrow(result_matrix), ncol=ncol(result_matrix))
correct_clustering <- as.vector(result_matrix)
# calc_hamming(unique(result_matrix, MARGIN = 1))


n <- length(unique(as.vector(result_matrix)))
if(n!=max(as.vector(result_matrix))){
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
rasterVis::levelplot(result_matrix, att = n,
                     # col.regions = rainbow(odd_number_larger),
                     colorkey = list(at = regions,
                                     # col=rainbow(odd_number_larger)[keep_these_colors],
                                     labels = list(at = middle_of_regions, labels = as.character(1:n))),
                     xlab = 'Cells',
                     ylab = 'Target genes',
                     main='Generated data')


###########################################
############ TODO plot data ###############
###########################################



######################################################################
############ Generate further complicated examples ###################
######################################################################




# Set seed for example
set.seed(1234)

# Set variables ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# preallocate data list
num_iter <- 20
generated_complicated_data_list <- vector(mode = "list", length = num_iter)



if (!file.exists(file.path(output_path, "env_sim_complicated_data_data_biclust.rds")) |
    redo_flag) {

  for(iter in 1:num_iter){

    n_cell_clusters <- 2
    n_target_gene_clusters <- c(4, 2)  # Number of target gene clusters in each cell cluster
    n_target_genes <- 100
    n_regulator_genes <- 10
    n_cells <- c(200, 200)
    regulator_means <- c(0, 0) # For generating dummy data, regulator mean in each cell cluster
    coefficient_means <- list(c(1, 3, 5, 7), c(10, 20))  # For generating dummy data, coefficient means in each cell cluster
    coefficient_sds <- list(c(0.01, 0.01, 0.01, 0.01), c(0.01, 0.01))
    disturbed_fraction <- 0.1  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
    testing_penalization_data_gen <- c(0.1, 0.5)



    generated_complicated_data <- generate_dummy_data_for_cell_clustering(
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
      plot_suffix = "Simple",
      testing_penalization = testing_penalization_data_gen,
      generate_counts             = FALSE,
      check_results               = FALSE,
      trivial_regulator_networks  = FALSE,
      pearson_regulators          = TRUE
    )

    generated_complicated_data_list[[iter]] <- generated_complicated_data

  }




  saveRDS(generated_complicated_data_list,
          file.path(output_path, "env_sim_complicated_data_data_biclust.rds"))

} else {
  generated_complicated_data_list <- readRDS(file.path(output_path, "env_sim_complicated_data_data_biclust.rds"))

}


