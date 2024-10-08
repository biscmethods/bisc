


if (!file.exists(file.path(output_path, "biclustbiclust_results_list.rds")) |
    redo_flag) {

  biclustbiclust_results_list <- vector(mode = "list", length = length(generated_complicated_data_list))

  for(iter in 1:length(generated_complicated_data_list )){

    print(paste0('Now running outer iteration ', iter))

    biclustbiclust_results_list[[iter]] <- bicluctbiclust(
      data = as.matrix(generated_complicated_data_list[[iter]]$dat[, 1:n_target_genes])
    )

  }

  saveRDS(
    biclustbiclust_results_list,
    file.path(output_path, "biclustbiclust_results_list.rds")
  )

} else {

  biclustbiclust_results_list    <- readRDS(file.path(output_path, "biclustbiclust_results_list.rds"))

}






