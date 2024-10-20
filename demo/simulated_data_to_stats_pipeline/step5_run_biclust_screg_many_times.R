
# Run model fits

if (!file.exists(file.path(output_path, "biclust_screg_results_list.rds")) |
    redo_flag) {


  biclust_screg_results_list <- vector(mode = "list", length = length(scenarios))

  for(iter in 1:length(biclust_screg_results_list )){

    cat("\n")
    print(paste0('Now running outer iteration ', iter, '/', length(scenarios), ' in step 5.'))
    cat("\n")
    set.seed(12)
    biclust_screg_results_list[[iter]] <- biclustscreg_iteration(plot_heatmap = FALSE,
                                                                 plot_title = paste0("heatmap_biclustscreg_", iter),
                                                                 penalization_lambdas = c(0.2), # c( 0.00001, 0.1, 0.2, 0.5)
                                                                 biclustscreg_results = NULL, # You can feed old results or calculate new ones
                                                                 cell_id = scenarios[[iter]]$cell_id,
                                                                 biclust_input_data = scenarios[[iter]]$biclust_input_data,
                                                                 output_path,  # Output path for biclust_screg for alluvial plots etc
                                                                 n_target_genes = scenarios[[iter]]$n_target_genes,
                                                                 n_total_cells = scenarios[[iter]]$n_total_cells,
                                                                 n_target_gene_clusters = scenarios[[iter]]$n_target_gene_clusters,
                                                                 n_cell_clusters = scenarios[[iter]]$n_cell_clusters,
                                                                 ind_targetgenes = scenarios[[iter]]$ind_targetgenes,
                                                                 ind_reggenes = scenarios[[iter]]$ind_reggenes,
                                                                 generated_data = scenarios[[iter]]$generated_data,
                                                                 correct_clustering = scenarios[[iter]]$correct_clustering,  # The correct biclustering (one unique number for each gene module)
                                                                 disturbed_initial_cell_clust = scenarios[[iter]]$disturbed_initial_cell_clust,
                                                                 itercap = 20)

  }

  saveRDS(
    biclust_screg_results_list,
    file.path(output_path, "biclust_screg_results_list.rds")
  )

} else {

  biclust_screg_results_list    <- readRDS(file.path(output_path, "biclust_screg_results_list.rds"))

}



