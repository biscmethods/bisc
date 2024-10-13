# !/usr/bin/Rscript

source(file.path(R_path, "biclust_scregclust.R"))


biclustscreg_iteration <- function(plot_heatmap=FALSE,
                                   penalization_lambdas = c(0.2), # c( 0.00001, 0.1, 0.2, 0.5)
                                   biclustscreg_results = NULL,  # You can feed old results or calculate new ones
                                   cell_id,
                                   biclust_input_data,
                                   output_path,  # Output path for biclust_screg for alluvial plots etc
                                   n_target_genes,
                                   n_total_cells,
                                   n_target_gene_clusters,
                                   n_cell_clusters,
                                   ind_targetgenes,
                                   ind_reggenes,
                                   generated_data,
                                   correct_clustering,  # The correct biclustering (one unique number for each gene module)
                                   disturbed_initial_cell_clust){
  # Run biclust_screg -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  if(is.null(biclustscreg_results)){
    biclustscreg_results <- vector(mode = "list", length = length(penalization_lambdas))

    for (i_penalization_lambda in seq_along(penalization_lambdas)) {
      print("", quote = FALSE)
      print(paste("Running biclust for penalization_lambda", penalization_lambdas[i_penalization_lambda]), quote = FALSE)

      biclustscreg_results[[i_penalization_lambda]] <- biclust_scregclust(
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
        plot_suffix                  = '',
        always_use_flat_prior        = FALSE,
        use_garbage_cluster_targets  = FALSE
      )
    }
  }

  print("", quote = FALSE)
  print("", quote = FALSE)
  print("Biclust results:", quote = FALSE)
  for (i_penalization_lambda in seq_along(penalization_lambdas)) {
    if (is.na(biclustscreg_results[i_penalization_lambda])) {
      print(paste("penalization_lambda", penalization_lambdas[i_penalization_lambda], "is NA"), quote = FALSE)
    } else if (is.null(biclustscreg_results[i_penalization_lambda])) {
      print(paste("penalization_lambda", penalization_lambdas[i_penalization_lambda], "is NULL"), quote = FALSE)
    } else {
      print(paste("penalization_lambda", penalization_lambdas[i_penalization_lambda], "is ok with rand index", biclustscreg_results[[i_penalization_lambda]]$rand_index), quote = FALSE)
    }
  }

  # Calculate RIs and maybe plot heatmaps ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  if(plot_heatmap){
    plots <- vector(mode = "list", length = length(penalization_lambdas))
  }
  RIs <- vector(mode = "list", length = length(penalization_lambdas))

  for(i_res in 1:length(penalization_lambdas)){
    current_biclust_result <- biclustscreg_results[[i_res]]
    if(is.list(current_biclust_result)){ # For no result current_biclust_result is just NA and then we can't calculate RI
      cell_cluster_allocation <- current_biclust_result$cell_cluster_allocation
      target_gene_allocation <- current_biclust_result$scregclust_final_result_module


      # Assign one unique number to each gene module --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
      biclust_result_matrix <- matrix(0, nrow = n_total_cells, ncol = n_target_genes)

      # Check each entry in the matrix (each cell and gene pair), and assign a string-number to each unique "cell-cluster-gene-module".
      for (i in 1:n_total_cells) {
        cluster <- cell_cluster_allocation[i]
        gene_allocation <- target_gene_allocation[[cluster]][1:n_target_genes]
        gene_allocation[gene_allocation==-1] <- 0
        # Create unique numbers for each pair
        biclust_result_matrix[i, ] <- paste0(cluster, gene_allocation)
      }

      # Convert the string-numbers to numeric matrix (starting from 1 this time)
      biclust_result_matrix <- matrix(as.numeric(biclust_result_matrix), nrow = n_total_cells, ncol = n_target_genes)
      biclust_result_matrix <- matrix(as.integer(as.factor(biclust_result_matrix)), nrow=nrow(biclust_result_matrix), ncol=ncol(biclust_result_matrix))


      # Calculate RIs ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
      true_cell_cluster_allocation <- generated_data$true_cell_clust
      true_target_gene_allocation <- generated_data$true_target_gene_allocation
      RI_cell_clustering_biclustscreg <- round(aricode::RI(cell_cluster_allocation, true_cell_cluster_allocation), 2)
      print(paste("Lambda", penalization_lambdas[i_res]), quote=FALSE)
      print(" Cell clustering RI for biclust_screg", quote=FALSE)
      print(paste("  ", RI_cell_clustering_biclustscreg), quote=FALSE)
      print(" Gene module clustering RI for biclust_screg", quote=FALSE)
      RI_gene_clustering_biclustscreg_all_string <- ""
      RI_gene_clustering_biclustscreg_all <- vector(length = n_cell_clusters)
      for(i_cell_cluster in 1:length(target_gene_allocation)){
        RI_gene_clustering_biclustscreg <- round(aricode::RI(target_gene_allocation[[i_cell_cluster]][1:n_target_genes], true_target_gene_allocation[[i_cell_cluster]][1:n_target_genes]), 2)
        print(paste("  For cell cluster", i_cell_cluster,":", RI_gene_clustering_biclustscreg), quote=FALSE)
        RI_gene_clustering_biclustscreg_all_string <- paste(RI_gene_clustering_biclustscreg_all_string, RI_gene_clustering_biclustscreg, sep=" ")
        RI_gene_clustering_biclustscreg_all[i_cell_cluster] <- RI_gene_clustering_biclustscreg
      }
      print(" Bicluster RI för biclust_screg",quote=FALSE)
      RI_biclust_biclustscreg <- round(aricode::RI(as.vector(biclust_result_matrix), correct_clustering), 2)
      print(paste("  ", RI_biclust_biclustscreg), quote=FALSE)

      # Save RIs
      RIs[[i_res]] <- list("penalization_lambda" = penalization_lambdas[i_res],
                           "RI_cell_clustering_biclustscreg" = RI_cell_clustering_biclustscreg,
                           "RI_gene_clustering_biclustscreg" = RI_gene_clustering_biclustscreg_all,
                           "RI_biclust_biclustscreg" = RI_biclust_biclustscreg)


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
                                               main=paste0('biclust_screg, lambda: ', penalization_lambdas[i_res],
                                                           "\nCell cluster RI:", RI_cell_clustering_biclustscreg,
                                                           "\nGene modules RI: ", RI_gene_clustering_biclustscreg_all_string,
                                                           "\nBiclust RI:", RI_biclust_biclustscreg))
      }

    }else{
      # Save RIs
      RIs[[i_res]] <- list("penalization_lambda" = penalization_lambdas[i_res],
                           "RI_cell_clustering_biclustscreg" = NA,
                           "RI_gene_clustering_biclustscreg" = NA,
                           "RI_biclust_biclustscreg" = NA)
    }
  }

  if(plot_heatmap){
    # Plot in IDE
    plot_construction <- cowplot::plot_grid(plotlist = plots,  align = 'vh', axis = 'tblr')
    print(plot_construction)
    # Plot to file
    png(file.path(output_path, paste0("biclustscreg_heatmap.png")), width = 1024, height = 480, units = "px")
    print(cowplot::plot_grid(plotlist = plots,  align = 'vh', axis = 'tblr'))
    dev.off()
  }

  return(list("RIs"=RIs,
              "biclustscreg_results"=biclustscreg_results)
  )
}



# Example use (run step1 first) -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Runs only when script is run by itself
# || interactive()
if (sys.nframe() == 0) {
  # Set seed for example
  set.seed(1234)
  res <- biclustscreg_iteration(plot_heatmap = FALSE,
                                penalization_lambdas = c(0.2, 1.0), # c( 0.00001, 0.1, 0.2, 0.5)
                                biclustscreg_results = NULL, # You can feed old results or calculate new ones
                                cell_id = scenarios[[1]]$cell_id,
                                biclust_input_data = scenarios[[1]]$biclust_input_data,
                                output_path,  # Output path for biclust_screg for alluvial plots etc
                                n_target_genes = scenarios[[1]]$n_target_genes,
                                n_total_cells = scenarios[[1]]$n_total_cells,
                                n_target_gene_clusters = scenarios[[1]]$n_target_gene_clusters,
                                n_cell_clusters = scenarios[[1]]$n_cell_clusters,
                                ind_targetgenes = scenarios[[1]]$ind_targetgenes,
                                ind_reggenes = scenarios[[1]]$ind_reggenes,
                                generated_data = scenarios[[1]]$generated_data,
                                correct_clustering = scenarios[[1]]$correct_clustering,  # The correct biclustering (one unique number for each gene module)
                                disturbed_initial_cell_clust = scenarios[[1]]$disturbed_initial_cell_clust)

}

