# !/usr/bin/Rscript

if (!file.exists(file.path(output_path, "bisc_results_list.rds")) | redo_flag) {

  bisc_results_list <- vector(mode = "list", length = length(scenarios))

  for(iter in 1:length(bisc_results_list)){
    cat("\n")
    print(paste0('Now running outer iteration ', iter, '/', length(scenarios), ' in step 5.'))
    cat("\n")
    bisc_results_list[[iter]] <- bisc_iteration(plot_heatmap = FALSE,
                                                plot_title = paste0("heatmap_bisc_", iter),
                                                penalization_lambdas = c(0.2), # c( 0.00001, 0.1, 0.2, 0.5)
                                                bisc_results = NULL, # You can feed old results or calculate new ones
                                                cell_id = scenarios[[iter]]$cell_id,
                                                biclust_input_data = scenarios[[iter]]$biclust_input_data,
                                                output_path,  # Output path for bisc for alluvial plots etc
                                                n_target_genes = scenarios[[iter]]$n_target_genes,
                                                n_total_cells = scenarios[[iter]]$n_total_cells,
                                                n_target_gene_clusters = scenarios[[iter]]$n_target_gene_clusters,
                                                n_cell_clusters = scenarios[[iter]]$n_cell_clusters,
                                                ind_targetgenes = scenarios[[iter]]$ind_targetgenes,
                                                ind_reggenes = scenarios[[iter]]$ind_reggenes,
                                                generated_data = scenarios[[iter]]$generated_data,
                                                correct_clustering = scenarios[[iter]]$correct_clustering,  # The correct biclustering (one unique number for each gene module)
                                                disturbed_initial_cell_clust = scenarios[[iter]]$disturbed_initial_cell_clust,
                                                itercap = 20,
                                                biclust_input_data_test = scenarios[[iter]]$biclust_input_data_test,
                                                n_total_cells_test = scenarios[[iter]]$n_total_cells_test,
                                                correct_clustering_test = scenarios[[iter]]$correct_clustering_tes,
                                                seeds = seq(10))

  }
  saveRDS(
    bisc_results_list,
    file.path(output_path, "bisc_results_list.rds")
  )
} else {
  bisc_results_list    <- readRDS(file.path(output_path, "bisc_results_list.rds"))
}



