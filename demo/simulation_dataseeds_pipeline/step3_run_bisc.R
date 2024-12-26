# !/usr/bin/Rscript

source(file.path(R_path, "bisc.R"))
source(file.path(R_path, "bisc_predict.R"))


bisc_iteration <- function(plot_heatmap=FALSE,
                           plot_title = "bisc_heatmap",
                           penalization_lambdas, # c( 0.00001, 0.1, 0.2, 0.5)
                           bisc_results = NULL,  # You can feed old results or calculate new ones
                           cell_id,
                           biclust_input_data,
                           output_path,  # Output path for bisc for alluvial plots etc
                           n_target_genes,
                           n_total_cells,
                           n_target_gene_clusters,
                           n_cell_clusters,
                           ind_targetgenes,
                           ind_reggenes,
                           generated_data,
                           correct_clustering,  # The correct biclustering (one unique number for each gene module)
                           disturbed_initial_cell_clust,
                           itercap,
                           biclust_input_data_test,
                           n_total_cells_test,
                           correct_clustering_test,
                           seeds){
  # Run bisc -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  if(is.null(bisc_results)){
    bisc_results <- vector(mode = "list", length = length(penalization_lambdas)*length(seeds))
    info <- vector(mode = "list", length = length(penalization_lambdas)*length(seeds))
    i <- 0
    for(current_seed in seeds){
      set.seed(current_seed)
      for (current_penalization_lambda in penalization_lambdas) {
        i <- i + 1
        print("", quote = FALSE)
        print(paste("Running biclust for penalization_lambda", current_penalization_lambda), quote = FALSE)
        plot_suffix <- paste0("_lambda_", current_penalization_lambda,
                              "_nCellClusters_", n_cell_clusters,
                              "_nRegulatorGenes_", length(ind_reggenes),
                              "_nTargetGenes_", n_target_genes,
                              "_nCells_", n_total_cells,
                              "_seed_", current_seed)
        info[[i]]$penalization_lambda <- current_penalization_lambda
        info[[i]]$seed <- current_seed
        bisc_results[[i]] <- bisc(
          dat = biclust_input_data,
          cell_id = cell_id,
          true_cell_cluster_allocation = factor(generated_data$true_cell_clust),
          max_iter = itercap,
          n_target_gene_clusters = n_target_gene_clusters,
          initial_clustering = disturbed_initial_cell_clust,
          n_cell_clusters = n_cell_clusters,
          ind_targetgenes = ind_targetgenes,
          ind_reggenes = ind_reggenes,
          output_path = output_path,
          penalization_lambda = current_penalization_lambda,
          calculate_optimization_target = TRUE,
          calculate_silhoutte = TRUE,
          calculate_davies_bouldin_index = TRUE,
          plot_suffix                  = plot_suffix,
          always_use_flat_prior        = FALSE,
          use_garbage_cluster_targets  = FALSE
        )
      }
    }
  }

  print("", quote = FALSE)
  print("", quote = FALSE)
  print(" Biclust results:", quote = FALSE)
  for (i in seq_along(info)) {
    if (is.na(bisc_results[i])) {
      print(paste("  penalization_lambda", info[[i]]$penalization_lambda, "is NA"), quote = FALSE)
    } else if (is.null(bisc_results[i])) {
      print(paste("  penalization_lambda",info[[i]]$penalization_lambda, "is NULL"), quote = FALSE)
    } else {
      print(paste("  penalization_lambda", info[[i]]$penalization_lambda, "is ok with rand index", bisc_results[[i]]$rand_index), quote = FALSE)
    }
  }

  # Calculate RIs and maybe plot heatmaps ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  if(plot_heatmap){
    plots <- vector(mode = "list", length = length(info))
  }
  RIs <- vector(mode = "list", length = length(info))

  for(i_res in seq_along(info)){
    current_biclust_result <- bisc_results[[i_res]]
    if(is.list(current_biclust_result)){ # For no result current_biclust_result is just NA and then we can't calculate RI

      cell_cluster_allocation_train <- current_biclust_result$cell_cluster_allocation
      cell_cluster_allocation <- bisc_predict(new_data = biclust_input_data_test,     # as in scenarios[[1]]$biclust_input_data
                                              fitted_model = bisc_results[[i_res]], # as in bisc_results_list[[1]]$bisc_results[[1]]
                                              prior_cluster_proportions = NULL,
                                              calculate_optimization_target = FALSE,
                                              seed = info[[i_res]]$penalization_lambda
      )$cell_cluster_allocation

      # logic check
      if(length(cell_cluster_allocation) != n_total_cells_test) {
        cat("\n")
        cat('Cluster allocation on test set of wrong length.')
        cat("\n")
      }


      target_gene_allocation <- current_biclust_result$scregclust_final_result_module


      # Assign one unique number to each gene module --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
      biclust_result_matrix <- matrix(0, nrow = n_total_cells_test, ncol = n_target_genes)

      # Check each entry in the matrix (each cell and gene pair),
      # and assign a string-number to each unique "cell-cluster-gene-module".
      for (i in 1:n_total_cells_test) {

        if(is.na(cell_cluster_allocation[i])){
          cat("Cell cluster allocation is NA. Random cluster will be assigned.\n")
          cell_cluster_allocation[i] <- sample.int(length(target_gene_allocation),1)
        }

        cluster <- cell_cluster_allocation[i]
        gene_allocation <- target_gene_allocation[[cluster]][1:n_target_genes]
        gene_allocation[gene_allocation==-1] <- 0
        # Create unique numbers for each pair
        biclust_result_matrix[i, ] <- paste0(cluster, gene_allocation)
      }

      # Convert the string-numbers to numeric matrix (starting from 1 this time)
      biclust_result_matrix <- matrix(as.numeric(biclust_result_matrix),
                                      nrow = n_total_cells_test,
                                      ncol = n_target_genes)

      biclust_result_matrix <- matrix(as.integer(as.factor(biclust_result_matrix)),
                                      nrow=nrow(biclust_result_matrix),
                                      ncol=ncol(biclust_result_matrix))


      # Calculate RIs ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
      true_cell_cluster_allocation_train <- generated_data$true_cell_clust
      true_cell_cluster_allocation <- generated_data$true_cell_clust_test

      correct_bicluster_test <- correct_clustering_test

      true_target_gene_allocation <- generated_data$true_target_gene_allocation #same for test and train data, no change


      RI_cell_clustering_bisc <- aricode::RI(cell_cluster_allocation, true_cell_cluster_allocation)
      print(paste(" Lambda", info[[i_res]]$penalization_lambda), quote=FALSE)
      print("  Cell clustering RI for bisc", quote=FALSE)
      print(paste("  ", RI_cell_clustering_bisc), quote=FALSE)
      print("  Gene module clustering RI for bisc", quote=FALSE)
      RI_gene_clustering_bisc_all_string <- ""
      RI_gene_clustering_bisc_all <- vector(length = n_cell_clusters)
      for(i_cell_cluster in seq_along(target_gene_allocation)){
        RI_gene_clustering_bisc <- aricode::RI(target_gene_allocation[[i_cell_cluster]][1:n_target_genes], true_target_gene_allocation[[i_cell_cluster]][1:n_target_genes])
        print(paste("   For cell cluster", i_cell_cluster,":", RI_gene_clustering_bisc), quote=FALSE)
        RI_gene_clustering_bisc_all_string <- paste(RI_gene_clustering_bisc_all_string, RI_gene_clustering_bisc, sep=" ")
        RI_gene_clustering_bisc_all[i_cell_cluster] <- RI_gene_clustering_bisc
      }
      print("  Bicluster RI för bisc",quote=FALSE)
      RI_biclust_bisc <- aricode::RI(as.vector(biclust_result_matrix), correct_bicluster_test)
      print(paste("  ", RI_biclust_bisc), quote=FALSE)

      # Save RIs
      RIs[[i_res]] <- list("seed" = info[[i_res]]$seed,
                           "penalization_lambda" = info[[i_res]]$penalization_lambda,
                           "RI_cell_clustering_bisc" = RI_cell_clustering_bisc,
                           "RI_gene_clustering_bisc" = RI_gene_clustering_bisc_all,
                           "RI_biclust_bisc" = RI_biclust_bisc)


      # Potentially construct plots -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
      if(plot_heatmap){
        n <- length(unique(as.vector(biclust_result_matrix)))
        regions <- seq(1, n, length.out = n + 1)
        middle_of_regions <- (regions[-1] + regions[-length(regions)]) / 2
        odd_number_larger <- ifelse(n %% 2 == 0, n + 1, n)
        if(n %% 2 == 1){
          keep_these_colors = 1:n
        }else{
          keep_these_colors <- setdiff(1:odd_number_larger, (odd_number_larger + 1) / 2 + 1)
        }
        plots[[i_res]] <- rasterVis::levelplot(biclust_result_matrix,
                                               att = n,
                                               # col.regions = rainbow(odd_number_larger),
                                               colorkey = list(at = regions,
                                                               # col=rainbow(odd_number_larger)[keep_these_colors],
                                                               labels = list(at = middle_of_regions, labels = as.character(1:n))),
                                               xlab = 'Cells',
                                               ylab = 'Target genes',
                                               main=paste0('bisc, lambda: ', info[[i_res]]$penalization_lambda,
                                                           "\nSeed:", info[[i_res]]$seed,
                                                           "\nCell cluster RI:", RI_cell_clustering_bisc,
                                                           "\nGene modules RI: ", RI_gene_clustering_bisc_all_string,
                                                           "\nBiclust RI:", RI_biclust_bisc,
                                                           "\nConverged:", bisc_results[[i_res]]$converged))
      }

    }else{
      # Save RIs
      RIs[[i_res]] <- list("seed" = info[[i_res]]$seed,
                           "penalization_lambda" = info[[i_res]]$penalization_lambda,
                           "RI_cell_clustering_bisc" = NA,
                           "RI_gene_clustering_bisc" = NA,
                           "RI_biclust_bisc" = NA)
    }
  }

  if(plot_heatmap){
    # Plot in IDE
    plot_construction <- cowplot::plot_grid(plotlist = plots,  align = 'vh', axis = 'tblr')
    print(plot_construction)
    # Plot to file
    png(file.path(output_path, paste0(plot_title, ".png")), width = 1024, height = 480, units = "px")
    print(cowplot::plot_grid(plotlist = plots,  align = 'vh', axis = 'tblr'))
    dev.off()
  }

  return(list("RIs"=RIs,
              "bisc_results"=bisc_results)
  )
}



# Example use (run step1 first) -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Runs only when script is run by itself
# || interactive()
if (sys.nframe() == 0) {
  # Set seed for example
  res <- bisc_iteration(plot_heatmap = TRUE,
                        plot_title = "heatmap_bisc_lambda_0.2",
                        penalization_lambdas = c(0.2), # c( 0.00001, 0.1, 0.2, 0.5)
                        bisc_results = NULL, # You can feed old results or calculate new ones
                        cell_id = scenarios[[1]]$cell_id,
                        biclust_input_data = scenarios[[1]]$biclust_input_data,
                        output_path,  # Output path for bisc for alluvial plots etc
                        n_target_genes = scenarios[[1]]$n_target_genes,
                        n_total_cells = scenarios[[1]]$n_total_cells,
                        n_target_gene_clusters = scenarios[[1]]$n_target_gene_clusters,
                        n_cell_clusters = scenarios[[1]]$n_cell_clusters,
                        ind_targetgenes = scenarios[[1]]$ind_targetgenes,
                        ind_reggenes = scenarios[[1]]$ind_reggenes,
                        generated_data = scenarios[[1]]$generated_data,
                        correct_clustering = scenarios[[1]]$correct_clustering,  # The correct biclustering (one unique number for each gene module)
                        disturbed_initial_cell_clust = scenarios[[1]]$disturbed_initial_cell_clust,
                        itercap=20,
                        biclust_input_data_test = scenarios[[1]]$biclust_input_data_test,
                        n_total_cells_test = scenarios[[1]]$n_total_cells_test,
                        correct_clustering_test = scenarios[[1]]$correct_clustering_test,
                        seeds=1234)

  res <- bisc_iteration(plot_heatmap = TRUE,
                        plot_title = "heatmap_bisc_lambda_1.0",
                        penalization_lambdas = c(1.0), # c( 0.00001, 0.1, 0.2, 0.5)
                        bisc_results = NULL, # You can feed old results or calculate new ones
                        cell_id = scenarios[[1]]$cell_id,
                        biclust_input_data = scenarios[[1]]$biclust_input_data,
                        output_path,  # Output path for bisc for alluvial plots etc
                        n_target_genes = scenarios[[1]]$n_target_genes,
                        n_total_cells = scenarios[[1]]$n_total_cells,
                        n_target_gene_clusters = scenarios[[1]]$n_target_gene_clusters,
                        n_cell_clusters = scenarios[[1]]$n_cell_clusters,
                        ind_targetgenes = scenarios[[1]]$ind_targetgenes,
                        ind_reggenes = scenarios[[1]]$ind_reggenes,
                        generated_data = scenarios[[1]]$generated_data,
                        correct_clustering = scenarios[[1]]$correct_clustering,  # The correct biclustering (one unique number for each gene module)
                        disturbed_initial_cell_clust = scenarios[[1]]$disturbed_initial_cell_clust,
                        itercap=20,
                        biclust_input_data_test = scenarios[[1]]$biclust_input_data_test,
                        n_total_cells_test = scenarios[[1]]$n_total_cells_test,
                        correct_clustering_test = scenarios[[1]]$correct_clustering_test,
                        seeds=1234)


  # Run for all scenarios
  # res <- vector(mode = "list", length = length(scenarios))
  # for(i in seq_along(scenarios)){
  #   set.seed(12)
  #   res[[i]] <- bisc_iteration(plot_heatmap = FALSE,
  #                                      plot_title = "bisc_heatmap",
  #                                      penalization_lambdas = c(0.2), # c( 0.00001, 0.1, 0.2, 0.5)
  #                                      bisc_results = NULL, # You can feed old results or calculate new ones
  #                                      cell_id = scenarios[[i]]$cell_id,
  #                                      biclust_input_data = scenarios[[i]]$biclust_input_data,
  #                                      output_path,  # Output path for bisc for alluvial plots etc
  #                                      n_target_genes = scenarios[[i]]$n_target_genes,
  #                                      n_total_cells = scenarios[[i]]$n_total_cells,
  #                                      n_target_gene_clusters = scenarios[[i]]$n_target_gene_clusters,
  #                                      n_cell_clusters = scenarios[[i]]$n_cell_clusters,
  #                                      ind_targetgenes = scenarios[[i]]$ind_targetgenes,
  #                                      ind_reggenes = scenarios[[i]]$ind_reggenes,
  #                                      generated_data = scenarios[[i]]$generated_data,
  #                                      correct_clustering = scenarios[[i]]$correct_clustering,  # The correct biclustering (one unique number for each gene module)
  #                                      disturbed_initial_cell_clust = scenarios[[i]]$disturbed_initial_cell_clust,
  #                                      itercap = 20)
  #
  # }
}
