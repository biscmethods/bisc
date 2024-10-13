# !/usr/bin/Rscript

library(biclust)

# Run biclust::biclust --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

bicluctbiclust <- function(data){ #centralise this function call so that we only need to set arguments once
  biclust::biclust(
    data,
    method=BCPlaid(),
    background=FALSE,
    iter.startup=100,
    iter.layer=100,
    back.fit=100,
    row.release=0.7,
    col.release=0.7,
    shuffle=100,
    max.layers=5,
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

  # Print the distance matrix
  return(dist_matrix)
}

cluster_biclust <- function(biclust_result, biclust_input_data, n_target_genes, n_total_cells, n_target_gene_clusters, n_cell_clusters){
  biclust_input_data <- as.matrix(biclust_input_data)[,1:n_target_genes]

  # Result of biclustering as a matrix
  biclust_results_matrix <- matrix(0, nrow = n_total_cells, ncol = n_target_genes)
  for (i_n in 1:biclust_result@Number) {
    a <- as.matrix((biclust_result@RowxNumber[, i_n, drop = FALSE] %*% biclust_result@NumberxCol[i_n, , drop = FALSE]) == 1)
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
  print("Cell clusters found in biclust", quote=FALSE)
  print(res_cell_cluster, quote=FALSE)

  res_gene_cluster_per_cell_cluster <- vector("list", length = max(unique_cell_clusters))
  # Gene clusters inside cell clusters  -------------------
  for(i_cell_cluster in 1:max(unique_cell_clusters)){
    ind_cell_cluster <- res_cell_cluster==i_cell_cluster
    unique_genes <- unique(biclust_results_matrix[ind_cell_cluster,], MARGIN = 2)
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
      inds <- which(apply(biclust_results_matrix[ind_cell_cluster,], 2, function(x) return(all(x == unique_genes[,i]))))
      res_gene_cluster[inds] <- unique_gene_clusters[i]
    }

    res_gene_cluster_per_cell_cluster[[i_cell_cluster]] <- res_gene_cluster
    # Print the cluster assignments
    print(paste("Gene modules found in biclust for cell cluster", i_cell_cluster), quote=FALSE)
    print(res_gene_cluster, quote=FALSE)
  }

  return(list("res_cell_cluster"       = res_cell_cluster,
              "res_gene_cluster"       = res_gene_cluster_per_cell_cluster,
              "biclust_results_matrix" = biclust_results_matrix))
}

get_stats_biclustbiclust <- function(biclust_result,
                                     biclust_input_data,
                                     n_target_genes,
                                     n_total_cells,
                                     n_target_gene_clusters,
                                     n_cell_clusters,
                                     generated_data,
                                     correct_clustering) {

  res_biclust_clustering <- cluster_biclust(biclust_result         = biclust_result,
                                            biclust_input_data     = biclust_input_data,
                                            n_target_genes         = n_target_genes,
                                            n_total_cells          = n_total_cells,
                                            n_target_gene_clusters = n_target_gene_clusters,
                                            n_cell_clusters        = n_cell_clusters)

  res_cell_cluster <- res_biclust_clustering$res_cell_cluster
  res_gene_cluster <- res_biclust_clustering$res_gene_cluster
  biclust_results_matrix <- res_biclust_clustering$biclust_results_matrix

  # Get RI for cell and gene clustering
  true_cell_cluster_allocation <- generated_data$true_cell_clust
  true_target_gene_allocation <- generated_data$true_target_gene_allocation
  RI_cell_clustering_biclustbiclust <- round(aricode::RI(res_cell_cluster, true_cell_cluster_allocation), 2)
  print("Cell clustering RI for biclust::biclust", quote=FALSE)
  print(paste(" ", RI_cell_clustering_biclustbiclust), quote=FALSE)
  print("Gene module clustering RI for biclust::biclust", quote=FALSE)
  RI_gene_clustering_biclustbiclust_all <- ""

  for(i_cell_cluster in 1:length(res_gene_cluster)){
    RI_gene_clustering_biclustbiclust <- round(aricode::RI(res_gene_cluster[[i_cell_cluster]], true_target_gene_allocation[[i_cell_cluster]][1:n_target_genes]), 2)
    print(paste(" For cell cluster", i_cell_cluster,":", RI_gene_clustering_biclustbiclust), quote=FALSE)
    RI_gene_clustering_biclustbiclust_all <- paste(RI_gene_clustering_biclustbiclust_all, RI_gene_clustering_biclustbiclust, sep=" ")
  }

  print("Bicluster RI för biclust::biclust",quote=FALSE)
  RI_biclust_biclustbiclust <- round(aricode::RI(as.vector(biclust_results_matrix), correct_clustering), 2)
  print(paste(" ", RI_biclust_biclustbiclust), quote=FALSE)

  return(list ("biclust_results_matrix" = biclust_results_matrix,
               "RI_cell_clustering_biclustbiclust" = RI_cell_clustering_biclustbiclust,
               "RI_gene_clustering_biclustbiclust_all" = RI_gene_clustering_biclustbiclust_all,
               "RI_biclust_biclustbiclust"  = RI_biclust_biclustbiclust)
  )

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


biclustbiclust_iteration <- function(biclust_input_data     = biclust_input_data,
                                     n_target_genes         = n_target_genes,
                                     n_total_cells          = n_total_cells,
                                     n_target_gene_clusters = n_target_gene_clusters,
                                     n_cell_clusters        = n_cell_clusters,
                                     generated_data         = generated_data,
                                     correct_clustering     = correct_clustering){
  res1 <- bicluctbiclust(
    data = as.matrix(biclust_input_data[, 1:n_target_genes])
  )

  stats_biclustbiclust <- get_stats_biclustbiclust(biclust_result         = res1,
                                                   biclust_input_data     = biclust_input_data,
                                                   n_target_genes         = n_target_genes,
                                                   n_total_cells          = n_total_cells,
                                                   n_target_gene_clusters = n_target_gene_clusters,
                                                   n_cell_clusters        = n_cell_clusters,
                                                   generated_data         = generated_data,
                                                   correct_clustering     = correct_clustering)

  return(stats_biclustbiclust)
}




# Example use (run step1 first) -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Runs only when script is run by itself
# || interactive()
if (sys.nframe() == 0) {
  # Set seed for example
  set.seed(1234)

  stats_biclustbiclust <- biclustbiclust_iteration(biclust_input_data     = scenarios[[1]]$biclust_input_data,
                                                   n_target_genes         = scenarios[[1]]$n_target_genes,
                                                   n_total_cells          = scenarios[[1]]$n_total_cells,
                                                   n_target_gene_clusters = scenarios[[1]]$n_target_gene_clusters,
                                                   n_cell_clusters        = scenarios[[1]]$n_cell_clusters,
                                                   generated_data         = scenarios[[1]]$generated_data,
                                                   correct_clustering     = scenarios[[1]]$correct_clustering)

  constructed_plots <- plot_biclust_heatmap(biclust_results_matrix                = stats_biclustbiclust$biclust_results_matrix,
                                            RI_cell_clustering_biclustbiclust     = stats_biclustbiclust$RI_cell_clustering_biclustbiclust,
                                            RI_gene_clustering_biclustbiclust_all = stats_biclustbiclust$RI_gene_clustering_biclustbiclust_all,
                                            RI_biclust_biclustbiclust             = stats_biclustbiclust$RI_biclust_biclustbiclust)

  print(constructed_plots)
  png(file.path(output_path, paste0("biclustbiclust_heatmap.png")),
      width = 1024, height = 480, units = "px")
  print(constructed_plots)
  dev.off()

}
