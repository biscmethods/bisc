# !/usr/bin/Rscript

# Run model fits
if (!file.exists(file.path(data_path, "biclustbiclust_results_list.rds")) | redo_flag) {

  biclustbiclust_results_list <- vector(mode = "list", length = length(scenarios))

  for(iter in seq_along(biclustbiclust_results_list)){

    cat("\n")
    print(paste0('Now running outer iteration ', iter, '/', length(scenarios), ' in step 4.'))
    cat("\n")

    biclustbiclust_results_list[[iter]] <- get_stats_biclustbiclust_seeds(biclust_input_data     = scenarios[[iter]]$biclust_input_data,
                                                                          n_target_genes         = scenarios[[iter]]$n_target_genes,
                                                                          ind_targetgenes        = scenarios[[iter]]$ind_targetgenes,
                                                                          n_total_cells          = scenarios[[iter]]$n_total_cells,
                                                                          n_target_gene_clusters = scenarios[[iter]]$n_target_gene_clusters,
                                                                          n_cell_clusters        = scenarios[[iter]]$n_cell_clusters,
                                                                          generated_data         = scenarios[[iter]]$generated_data,
                                                                          correct_clustering     = scenarios[[iter]]$correct_clustering,
                                                                          seeds = c(12),
                                                                          do_biclust_with_regulators = TRUE)

    # constructed_plots <- plot_biclust_heatmap(biclust_results_matrix                = stats_biclustbiclust$biclust_results_matrix,
    #                                           RI_cell_clustering_biclustbiclust     = stats_biclustbiclust$RI_cell_clustering_biclustbiclust,
    #                                           RI_gene_clustering_biclustbiclust_all = stats_biclustbiclust$RI_gene_clustering_biclustbiclust_all,
    #                                           RI_biclust_biclustbiclust             = stats_biclustbiclust$RI_biclust_biclustbiclust)
    #
    # print(constructed_plots)
    # png(file.path(output_path, paste0("heatmap_biclustbiclust_",iter,".png")),
    #     width = 1024, height = 480, units = "px")
    # print(constructed_plots)
    # dev.off()

  }

  saveRDS(biclustbiclust_results_list, file.path(data_path, "biclustbiclust_results_list.rds"))

} else {

  biclustbiclust_results_list    <- readRDS(file.path(data_path, "biclustbiclust_results_list.rds"))

}


