

biclustbiclust_results_list <- vector(mode = "list", length = length(generated_complicated_data_list))


for(iter in 1:length(generated_complicated_data_list )){

  biclustbiclust_results_list[[iter]] <- bicluctbiclust(
    data = as.matrix(generated_complicated_data_list[[iter]]$dat[, 1:n_target_genes])
  )

}





