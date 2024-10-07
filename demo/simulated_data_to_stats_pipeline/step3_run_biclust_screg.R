#!/usr/bin/Rscript
# rm(list = ls())

library(here)  # To work with paths
library(patchwork)
library(rasterVis)
library(cowplot)
sink()


# options(warn=2)  # To convert warning messages into error messages which display row of error. For debugging.

# Get absolute path where script is located, by using relative paths.
demo_path <- here::here("demo")
R_path <- here::here("R")
output_path <- here::here("demo/simulated_data_to_stats_pipeline")
path_data <- here::here('data')

source(file.path(R_path, "biclust_scregclust.R"))

# Set seed for example
set.seed(1234)



# Run biclust_screg -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# TODO: Which lambda result shall we use? All? The best knowing the solution? Best would be if we can get some selection criteria working but shall we skip that for now?

penalization_lambdas <- c( 0.1, 0.3, 0.5) # c( 0.00001, 0.1, 0.2, 0.5)
BICLUST_RESULTS <- vector(mode = "list", length = length(penalization_lambdas))

if (!file.exists(file.path(path_data, "env_sim_simple_nogarb_res_biclust_sc.rds")) |
    redo_flag) {
  for (i_penalization_lambda in seq_along(penalization_lambdas)) {
    print("", quote = FALSE)
    print(paste(
      "Running biclust for penalization_lambda",
      penalization_lambdas[i_penalization_lambda]
    ),
    quote = FALSE)
    BICLUST_RESULTS[[i_penalization_lambda]] <- biclust_scregclust(
      dat = biclust_input_data,
      cell_id = cell_id,
      true_cell_cluster_allocation = factor(generated_data$true_cell_clust),
      max_iter = 100,
      n_target_gene_clusters = n_target_gene_clusters,
      initial_clustering = disturbed_initial_cell_clust,
      n_cell_clusters = n_cell_clusters,
      ind_targetgenes = ind_targetgenes,
      ind_reggenes = ind_reggenes,
      output_path = output_path,
      penalization_lambda = penalization_lambdas[i_penalization_lambda],
      use_complex_cluster_allocation = FALSE,
      calculate_BIC = FALSE,
      calculate_silhoutte = FALSE,
      calculate_davies_bouldin_index = FALSE,
      plot_suffix                  = "step3_biclust_screg",
      always_use_flat_prior        = FALSE,
      use_garbage_cluster_targets  = FALSE
    )
  }

  saveRDS(
    BICLUST_RESULTS,
    file.path(path_data, "env_sim_simple_nogarb_res_biclust_sc.rds")
  )

} else {
  BICLUST_RESULTS <- readRDS(file.path(path_data, "env_sim_simple_nogarb_res_biclust_sc.rds"))

}

print("", quote = FALSE)
print("", quote = FALSE)
for (i_penalization_lambda in seq_along(penalization_lambdas)) {
  if (is.na(BICLUST_RESULTS[i_penalization_lambda])) {
    print(paste("penalization_lambda", penalization_lambdas[i_penalization_lambda], "is NA"),
          quote = FALSE)
  } else if (is.null(BICLUST_RESULTS[i_penalization_lambda])) {
    print(paste("penalization_lambda", penalization_lambdas[i_penalization_lambda], "is NULL"),
          quote = FALSE)
  } else {
    print(
      paste(
        "penalization_lambda",
        penalization_lambdas[i_penalization_lambda],
        "is ok with rand index",
        BICLUST_RESULTS[[i_penalization_lambda]]$rand_index
      ),
      quote = FALSE
    )
  }
}


# Construct heatmap for our biclust_screg -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


plots <- vector(mode = "list", length = length(penalization_lambdas))
for(i_res in 1:length(penalization_lambdas)){
  res <- BICLUST_RESULTS[[i_res]]
  if(is.list(res)){
    cell_cluster_allocation <- res$cell_cluster_allocation
    target_gene_allocation <- res$scregclust_final_result_module

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

    # CALC RIs
    true_cell_cluster_allocation <- generated_data$true_cell_clust
    true_target_gene_allocation <- generated_data$true_target_gene_allocation
    RI_cell_clustering_biclustscreg <- round(aricode::RI(cell_cluster_allocation, true_cell_cluster_allocation), 2)
    print(paste("Lambda", penalization_lambdas[i_res]), quote=FALSE)
    print(" Cell clustering RI for biclust_screg", quote=FALSE)
    print(paste("  ", RI_cell_clustering_biclustscreg), quote=FALSE)
    print(" Gene module clustering RI for biclust_screg", quote=FALSE)
    RI_gene_clustering_biclustscreg_all <- ""
    for(i_cell_cluster in 1:length(res_gene_cluster)){
      RI_gene_clustering_biclustscreg <- round(aricode::RI(target_gene_allocation[[i_cell_cluster]][1:n_target_genes], true_target_gene_allocation[[i_cell_cluster]][1:n_target_genes]), 2)
      print(paste("  For cell cluster", i_cell_cluster,":", RI_gene_clustering_biclustscreg), quote=FALSE)
      RI_gene_clustering_biclustscreg_all <- paste(RI_gene_clustering_biclustscreg_all, RI_gene_clustering_biclustscreg, sep=" ")
    }
    print(" Bicluster RI för biclust_screg",quote=FALSE)
    RI_biclust_biclustscreg <- round(aricode::RI(as.vector(result_matrix), correct_clustering), 2)
    print(paste("  ", RI_biclust_biclustscreg), quote=FALSE)



    n <- length(unique(as.vector(result_matrix)))
    regions <- seq(1, n, length.out = n + 1)
    middle_of_regions <- (regions[-1] + regions[-length(regions)]) / 2
    odd_number_larger <- ifelse(n %% 2 == 0, n + 1, n)
    if(n %% 2 == 1){
      keep_these_colors = 1:n
    }else{
      keep_these_colors <- setdiff(1:odd_number_larger, (odd_number_larger + 1) / 2 + 1)
    }
    plots[[i_res]] <- rasterVis::levelplot(result_matrix,
                                           att = n,
                                           # col.regions = rainbow(odd_number_larger),
                                           colorkey = list(at = regions,
                                                           # col=rainbow(odd_number_larger)[keep_these_colors],
                                                           labels = list(at = middle_of_regions, labels = as.character(1:n))),
                                           xlab = 'Cells',
                                           ylab = 'Target genes',
                                           main=paste0('biclust_screg, lambda: ', penalization_lambdas[i_res],
                                                       "\nCell cluster RI:", RI_cell_clustering_biclustscreg,
                                                       "\nGene modules RI: ", RI_gene_clustering_biclustscreg_all,
                                                       "\nBiclust RI:", RI_biclust_biclustscreg))
  }
}

cowplot::plot_grid(plotlist = plots,  align = 'vh', axis = 'tblr')


