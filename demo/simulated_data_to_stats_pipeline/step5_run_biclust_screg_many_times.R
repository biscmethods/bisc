





if (!file.exists(file.path(output_path, "biclust_screg_results_list.rds")) |
    redo_flag) {


  biclust_screg_results_list <- vector(mode = "list", length = length(generated_complicated_data_list))

  for(iter in 1:length(generated_complicated_data_list )){

    biclust_screg_results_list[[iter]] <- run_biclust_scregclust(
      dat = generated_complicated_data_list[[iter]]$dat
    )

  }

  saveRDS(
    biclust_screg_results_list,
    file.path(output_path, "biclust_screg_results_list.rds")
  )

} else {

  biclust_screg_results_list    <- readRDS(file.path(output_path, "biclust_screg_results_list.rds"))

}





