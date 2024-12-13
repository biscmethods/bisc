# !/usr/bin/Rscript

library(biclust)

# Run biclust::biclust --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

biclustbiclust <- function(data){ #centralise this function call so that we only need to set arguments once
  biclust::biclust(
    data,
    method=BCPlaid(),
    background=FALSE,
    verbose=FALSE)
}




# Calculate stats and plot ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


calc_hamming <- function(matrix_data, between_rows=TRUE){
  hamming_dist <- function(vec1, vec2) {
    sum(vec1 != vec2)
  }

  # Number of rows (vectors)
  if(between_rows){
    n <- nrow(matrix_data)
  }else{
    n <- ncol(matrix_data)
  }


  # Initialize a distance matrix
  dist_matrix <- matrix(0, n, n)

  if(n==1){
    dist_matrix[1,1] <- 0
  }else{

    if(between_rows){
      # Calculate pairwise Hamming distances between ROWS
      for (i in 1:(n-1)) {
        for (j in (i+1):n) {
          dist_matrix[i, j] <- hamming_dist(matrix_data[i, ], matrix_data[j, ])
          dist_matrix[j, i] <- dist_matrix[i, j]  # Symmetric matrix
        }
      }
    }else{
      # Calculate pairwise Hamming distances between COLS
      for (i in 1:(n-1)) {
        for (j in (i+1):n) {
          dist_matrix[i, j] <- hamming_dist(matrix_data[, i], matrix_data[, j])
          dist_matrix[j, i] <- dist_matrix[i, j]  # Symmetric matrix
        }
      }
    }
  }
  # Print the distance matrix
  return(dist_matrix)
}


# Since biclustbiclust output row cluster assingments and col cluster assignment
# this predict function does the same.
# It simply correlates the train and test data per row/col and takes
# the cluster assignment from the best correlating row/col.
# Fo us, rows are cells, and cols are genes.
# > str(biclust_result@RowxNumber)
# logi [1:400, 1:5] FALSE FALSE TRUE FALSE FALSE FALSE ...
# > str(biclust_result@NumberxCol)
# logi [1:5, 1:106] TRUE FALSE FALSE FALSE FALSE TRUE ...
# predict_biclustbiclust <- function(biclust_result, train_data, test_data){
#   if (!identical(dim(train_data), dim(test_data))) {
#     print("train_data and test_data needs to have exactly the same dimensions.")
#   }
#
#   row_correlations <- apply(train_data, 1, function(row1) {
#     apply(test_data, 1, function(row2) {
#       cor(row1, row2)
#     })
#   })
#
#   col_correlations <- apply(train_data, 2, function(col1) {
#     apply(test_data, 2, function(col2) {
#       cor(col1, col2)
#     })
#   })
#
#   # Get the index of the highest correlation for each row/col
#   highest_corr_row_indices <- apply(row_correlations, 1, which.max)
#   highest_corr_col_indices <- apply(col_correlations, 2, which.max)
#
#   return(list("row_cluster_assignments" = biclust_result@RowxNumber[highest_corr_row_indices,],
#               "col_cluster_assignments" = biclust_result@NumberxCol[,highest_corr_col_indices]))
# }
#
#
# predict_biclustbiclust(biclust_result = stats_biclustbiclust$biclust_results_matrix,
#                        train_data = biclust_input_data,
#                        test_data = biclust_input_data)

get_stats_biclustbiclust <- function(biclust_input_data,
                                     n_target_genes,
                                     ind_targetgenes,
                                     n_total_cells,
                                     n_target_gene_clusters,
                                     n_cell_clusters,
                                     generated_data,
                                     correct_clustering,
                                     seed,
                                     do_biclust_with_regulators = TRUE,
                                     include_regulators_in_results = FALSE) {
  set.seed(seed)
  org_ind_targetgenes <- ind_targetgenes
  if(include_regulators_in_results && do_biclust_with_regulators){
    org_n_target_genes <- n_target_genes
    ind_targetgenes <- 1:ncol(biclust_input_data)
    n_target_genes <- ncol(biclust_input_data)
  }


  if(do_biclust_with_regulators){
    biclust_result <- biclustbiclust(data = as.matrix(biclust_input_data))
  }else{
    biclust_result <- biclustbiclust(data = as.matrix(biclust_input_data[ind_targetgenes]))
  }


  # Cluster the results from biclust --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


  # Result of biclustering as a matrix
  biclust_results_matrix <- matrix(0, nrow = n_total_cells, ncol = n_target_genes)
  for (i_n in 1:biclust_result@Number) {
    a <- as.matrix((biclust_result@RowxNumber[, i_n, drop = FALSE] %*% biclust_result@NumberxCol[i_n, ind_targetgenes, drop = FALSE]) == 1)
    biclust_results_matrix[a] <- i_n
  }


  # Cell clustering  -------------------
  # E.g. 400 cells x 100 genes
  # MARGIN=1 keeps the second/column/gene dimension intact - meaning it's cell clustering
  # e.g. 8 cells with 100 genes each
  unique_cells <- unique(biclust_results_matrix, MARGIN = 1)
  dist_matrix <- calc_hamming(unique_cells, between_rows=TRUE)
  distance_object <- stats::as.dist(dist_matrix)

  # Perform hierarchical clustering
  hc <- stats::hclust(distance_object)

  # Plot the dendrogram to visualize the clustering
  # plot(hc, labels = rownames(unique_cells))
  unique_cell_clusters <- stats::cutree(hc, k = n_cell_clusters)
  # unique_cell_clusters[unique_cell_clusters>2] = 2

  res_cell_cluster <- vector(length=nrow(biclust_results_matrix))
  for (i in 1:nrow(unique_cells)) {
    inds <- which(apply(biclust_results_matrix, 1, function(x) return(all(x == unique_cells[i,]))))
    res_cell_cluster[inds] <- unique_cell_clusters[i]
  }

  # Print the cluster assignments
  print("Cell clusters found in biclust::biclust", quote=FALSE)
  print(res_cell_cluster, quote=FALSE)

  res_gene_cluster_per_cell_cluster <- vector("list", length = max(unique_cell_clusters))
  # Gene clusters inside cell clusters  -------------------
  for(i_cell_cluster in 1:max(unique_cell_clusters)){
    ind_cell_cluster <- res_cell_cluster==i_cell_cluster
    unique_genes <- unique(biclust_results_matrix[ind_cell_cluster, ,drop = FALSE], MARGIN = 2, drop = FALSE)
    if(length(unique(unique_genes))==1){
      res_gene_cluster <- rep(1, ncol(biclust_results_matrix))
    }else{
      dist_matrix <- calc_hamming(unique_genes, between_rows=FALSE)
      distance_object <- stats::as.dist(dist_matrix)

      # Perform hierarchical clustering
      hc <- stats::hclust(distance_object)

      # Plot the dendrogram to visualize the clustering
      # plot(hc, labels = rownames(unique_cells))
      unique_gene_clusters <- stats::cutree(hc, k = min(n_target_gene_clusters[i_cell_cluster], length(hc$order)))
      # unique_gene_clusters[unique_gene_clusters>2] = 4

      res_gene_cluster <- vector(length=ncol(biclust_results_matrix))
      for (i in 1:ncol(unique_genes)) {
        inds <- which(apply(biclust_results_matrix[ind_cell_cluster, , drop=FALSE], 2, function(x) return(all(x == unique_genes[,i]))))
        res_gene_cluster[inds] <- unique_gene_clusters[i]
      }
    }
    res_gene_cluster_per_cell_cluster[[i_cell_cluster]] <- res_gene_cluster
    # Print the cluster assignments
    print(paste("Gene modules found in biclust::biclust for cell cluster", i_cell_cluster), quote=FALSE)
    print(res_gene_cluster, quote=FALSE)
  }

  # Get RI for cell and gene clustering -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  true_cell_cluster_allocation <- generated_data$true_cell_clust
  true_target_gene_allocation <- generated_data$true_target_gene_allocation
  RI_cell_clustering_biclustbiclust <- aricode::RI(res_cell_cluster, true_cell_cluster_allocation)
  print("Cell clustering RI for biclust::biclust", quote=FALSE)
  print(paste(" ", RI_cell_clustering_biclustbiclust), quote=FALSE)
  print("Gene module clustering RI for biclust::biclust", quote=FALSE)
  RI_gene_clustering_biclustbiclust_all <- ""

  for(i_cell_cluster in seq_along(res_gene_cluster_per_cell_cluster)){
    RI_gene_clustering_biclustbiclust <- aricode::RI(res_gene_cluster_per_cell_cluster[[i_cell_cluster]][org_ind_targetgenes], true_target_gene_allocation[[i_cell_cluster]][org_ind_targetgenes])
    print(paste(" For cell cluster", i_cell_cluster,":", RI_gene_clustering_biclustbiclust), quote=FALSE)
    RI_gene_clustering_biclustbiclust_all <- paste(RI_gene_clustering_biclustbiclust_all, RI_gene_clustering_biclustbiclust, sep=" ")
  }

  print("Bicluster RI för biclust::biclust",quote=FALSE)
  RI_biclust_biclustbiclust <- aricode::RI(as.vector(biclust_results_matrix[,org_ind_targetgenes]), correct_clustering)
  print(paste(" ", RI_biclust_biclustbiclust), quote=FALSE)

  return(list ("biclust_results_matrix" = biclust_results_matrix,
               "res_gene_cluster_per_cell_cluster" = res_gene_cluster_per_cell_cluster,
               "res_cell_cluster" = res_cell_cluster,
               "RI_cell_clustering_biclustbiclust" = RI_cell_clustering_biclustbiclust,
               "RI_gene_clustering_biclustbiclust_all" = RI_gene_clustering_biclustbiclust_all,
               "RI_biclust_biclustbiclust"  = RI_biclust_biclustbiclust,
               "seed" = seed)
  )

}

get_stats_biclustbiclust_seeds <- function(biclust_input_data,
                                           n_target_genes,
                                           ind_targetgenes,
                                           n_total_cells,
                                           n_target_gene_clusters,
                                           n_cell_clusters,
                                           generated_data,
                                           correct_clustering,
                                           seeds,
                                           do_biclust_with_regulators = TRUE,
                                           include_regulators_in_results = FALSE) {
  res <- vector(mode = "list", length = length(seeds))
  for(i_seed in seq_along(seeds)){
    cat("\n")
    print(paste0(' Now running n cell clusters ', n_cell_clusters,' for outer seed ', i_seed, '/', length(seeds), ' in step 4.1.'))
    cat("\n")
    res[[i_seed]] <- get_stats_biclustbiclust(biclust_input_data=biclust_input_data,
                                              n_target_genes=n_target_genes,
                                              ind_targetgenes=ind_targetgenes,
                                              n_total_cells=n_total_cells,
                                              n_target_gene_clusters=n_target_gene_clusters,
                                              n_cell_clusters=n_cell_clusters,
                                              generated_data=generated_data,
                                              correct_clustering=correct_clustering,
                                              seed=seeds[i_seed],
                                              do_biclust_with_regulators = do_biclust_with_regulators,
                                              include_regulators_in_results = include_regulators_in_results)
  }
  return(res)

}

# Plot the biclustering results
# b <- raster::ratify(raster::raster(b))
plot_biclust_heatmap <- function(biclust_results_matrix,
                                 RI_cell_clustering_biclustbiclust,
                                 RI_gene_clustering_biclustbiclust_all,
                                 RI_biclust_biclustbiclust){

  # This is just to fix colors with a unique legend
  n_unique_biclusters <- length(unique(as.vector(biclust_results_matrix)))
  regions <- seq(1, n_unique_biclusters, length.out = n_unique_biclusters + 1)
  middle_of_regions <- (regions[-1] + regions[-length(regions)]) / 2
  odd_number_larger <- ifelse(n_unique_biclusters %% 2 == 0, n_unique_biclusters + 1, n_unique_biclusters)
  if(n_unique_biclusters %% 2 == 1){
    keep_these_colors = 1:n_unique_biclusters
  }else{
    keep_these_colors <- setdiff(1:odd_number_larger, (odd_number_larger + 1) / 2 + 1)
  }
  constructed_plots <- rasterVis::levelplot(biclust_results_matrix, att = n_unique_biclusters,
                                            colorkey = list(at = regions,
                                                            labels = list(at = middle_of_regions, labels = as.character(1:n_unique_biclusters))),
                                            xlab = 'Cells',
                                            ylab = 'Target genes',
                                            main=paste0('biclust::biclust.\nCell cluster RI:',RI_cell_clustering_biclustbiclust,
                                                        "\nGene modules RI:", RI_gene_clustering_biclustbiclust_all,
                                                        "\nBiclust RI:", RI_biclust_biclustbiclust))
  return(constructed_plots)
}






# Example use (run step1 first) -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Runs only when script is run by itself
# || interactive()
if (sys.nframe() == 0) {
  # Set seed for example
  set.seed(1234)

  stats_biclustbiclust <- get_stats_biclustbiclust(biclust_input_data     = scenarios[[10]]$biclust_input_data,
                                                   n_target_genes         = scenarios[[10]]$n_target_genes,
                                                   ind_targetgenes        = scenarios[[10]]$ind_targetgenes,
                                                   n_total_cells          = scenarios[[10]]$n_total_cells,
                                                   n_target_gene_clusters = scenarios[[10]]$n_target_gene_clusters,
                                                   n_cell_clusters        = scenarios[[10]]$n_cell_clusters,
                                                   generated_data         = scenarios[[10]]$generated_data,
                                                   correct_clustering     = scenarios[[10]]$correct_clustering,
                                                   do_biclust_with_regulators = TRUE,
                                                   include_regulators_in_results = FALSE)

  constructed_plots <- plot_biclust_heatmap(biclust_results_matrix                = stats_biclustbiclust$biclust_results_matrix,
                                            RI_cell_clustering_biclustbiclust     = stats_biclustbiclust$RI_cell_clustering_biclustbiclust,
                                            RI_gene_clustering_biclustbiclust_all = stats_biclustbiclust$RI_gene_clustering_biclustbiclust_all,
                                            RI_biclust_biclustbiclust             = stats_biclustbiclust$RI_biclust_biclustbiclust)

  print(constructed_plots)
  png(file.path(output_path, paste0("heatmap_biclustbiclust.png")), width = 1024, height = 480, units = "px")
  print(constructed_plots)
  dev.off()

}
