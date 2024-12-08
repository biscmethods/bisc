# !/usr/bin/Rscript
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Assumes you have run step 4 first
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

library(biclust)

# Run biclust::biclust --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

biclustbiclust <- function(data){ #centralise this function call so that we only need to set arguments once
  biclust::biclust(
    data,
    method           = BCPlaid(),
    cluster          ="b",
    fit.model        = y ~ m + a + b,
    background       = TRUE,
    background.layer = NA,
    background.df    = 1,
    # row.release      = 0.1,
    # col.release      = 0.2,
    row.release      = 0.3,
    col.release      = 0.3,
    shuffle          = 2,
    back.fit         = 2,
    max.layers       = 50,
    iter.startup     = 5,
    iter.layer       = 100,
    verbose          = FALSE
  )
}

#
#


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


# arguments for dev
# first_converged <- which(converged)[1]
# str(all_res[[first_converged]][[1]]$call)
# biclust_input_data = d
# n_target_genes     = length(all_res[[first_converged]][[1]]$call$ind_targetgenes)
# ind_targetgenes    = all_res[[first_converged]][[1]]$call$ind_targetgenes
# n_total_cells      = length(all_res[[first_converged]][[1]]$call$cell_id)
# n_target_gene_clusters = all_res[[first_converged]][[1]]$call$n_target_gene_clusters
# n_cell_clusters        = all_res[[first_converged]][[1]]$call$n_cell_clusters
# true_cell_cluster_allocation = all_res[[first_converged]][[1]]$call$true_cell_cluster_allocation
# # generated_data         = all_res[[first_converged]][[1]]$call
# # correct_clustering
# seed  = 1
# do_biclust_with_regulators = TRUE
# include_regulators_in_results = FALSE


get_stats_biclustbiclust <- function(biclust_input_data = d,
                                     n_target_genes,
                                     ind_targetgenes,
                                     n_total_cells,
                                     n_target_gene_clusters,
                                     n_cell_clusters,
                                     true_cell_cluster_allocation, #<- generated_data$true_cell_clust
                                     # true_target_gene_allocation, #  <- generated_data$true_target_gene_allocation
                                     # correct_clustering,
                                     seed = 1,
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
    biclust_result <- biclustbiclust(data = as.matrix(biclust_input_data[ind_targetgenes,]))
  }


  # Cluster the results from biclust --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


  # Result of biclustering as a matrix
  biclust_results_matrix <- matrix(0, nrow = n_total_cells, ncol = n_target_genes)
  biclust_results_raw_cell_clustring <- matrix(0, nrow = n_total_cells, ncol = biclust_result@Number)
  for (i_n in 1:biclust_result@Number) {
    a <- as.matrix((biclust_result@RowxNumber[, i_n, drop = FALSE] %*% biclust_result@NumberxCol[i_n, ind_targetgenes, drop = FALSE]) == 1)
    biclust_results_matrix[a] <- i_n
    # biclust_results_raw_cell_clustring[a] <- i_n
    # print(length(ind_targetgenes))
    # print(dim(a))
    # print(length(biclust_result@RowxNumber[, i_n, drop = FALSE]))
    # print(length(biclust_result@NumberxCol[i_n, , drop = FALSE]))
    # cat("\n\n")
  }

  #
  # print(biclust_result@Number)
  # print(biclust_result@RowxNumber[, i_n, drop = FALSE])
  # exit()
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
  # unique_cell_clusters[unique_cell_clusters>n_cell_clusters] = n_cell_clusters

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
    print(paste("Gene modules found in biclust::biclust for cell cluster", i_cell_cluster), quote=FALSE)
    print(res_gene_cluster, quote=FALSE)
  }

  # Get RI for cell and gene clustering -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # true_cell_cluster_allocation <- generated_data$true_cell_clust
  # true_target_gene_allocation <- generated_data$true_target_gene_allocation
  #

  RI_cell_clustering_biclustbiclust <- aricode::RI(res_cell_cluster, true_cell_cluster_allocation)
  print("Cell clustering RI for biclust::biclust", quote=FALSE)
  print(paste(" ", RI_cell_clustering_biclustbiclust), quote=FALSE)
  print("Gene module clustering RI for biclust::biclust", quote=FALSE)
  RI_gene_clustering_biclustbiclust_all <- ""

  # for(i_cell_cluster in seq_along(res_gene_cluster_per_cell_cluster)){
  #   RI_gene_clustering_biclustbiclust <- round(aricode::RI(res_gene_cluster_per_cell_cluster[[i_cell_cluster]][org_ind_targetgenes], true_target_gene_allocation[[i_cell_cluster]][org_ind_targetgenes]), 2)
  #   print(paste(" For cell cluster", i_cell_cluster,":", RI_gene_clustering_biclustbiclust), quote=FALSE)
  #   RI_gene_clustering_biclustbiclust_all <- paste(RI_gene_clustering_biclustbiclust_all, RI_gene_clustering_biclustbiclust, sep=" ")
  # }

  # print("Bicluster RI för biclust::biclust",quote=FALSE)
  # RI_biclust_biclustbiclust <- round(aricode::RI(as.vector(biclust_results_matrix[,org_ind_targetgenes]), correct_clustering), 2)
  # print(paste(" ", RI_biclust_biclustbiclust), quote=FALSE)

  return(list ("biclust_results_matrix" = biclust_results_matrix,
               "res_gene_cluster_per_cell_cluster" = res_gene_cluster_per_cell_cluster,
               "res_cell_cluster" = res_cell_cluster,
               "RI_cell_clustering_biclustbiclust" = RI_cell_clustering_biclustbiclust,
               # "RI_gene_clustering_biclustbiclust_all" = RI_gene_clustering_biclustbiclust_all,
               # "RI_biclust_biclustbiclust"  = RI_biclust_biclustbiclust,
               "seed" = seed)
  )

}

get_stats_biclustbiclust_seeds <- function(biclust_input_data,
                                           n_target_genes,
                                           ind_targetgenes,
                                           n_total_cells,
                                           n_target_gene_clusters,
                                           n_cell_clusters,
                                           true_cell_cluster_allocation, #<- generated_data$true_cell_clust
                                           # true_target_gene_allocation, #  <- generated_data$true_target_gene_allocation
                                           # correct_clustering,
                                           seeds,
                                           do_biclust_with_regulators = TRUE,
                                           include_regulators_in_results = FALSE) {

  res <- vector(mode = "list", length = length(seeds))

  for(i_seed in seq_along(seeds)){

    print(paste0("RUNNING BCPLAID ITERATION ", i_seed, " OF ", length(seeds), "."))

    res[[i_seed]] <-tryCatch(
      expr =  get_stats_biclustbiclust(biclust_input_data=biclust_input_data,
                                       n_target_genes=n_target_genes,
                                       ind_targetgenes=ind_targetgenes,
                                       n_total_cells=n_total_cells,
                                       n_target_gene_clusters=n_target_gene_clusters,
                                       n_cell_clusters=n_cell_clusters,
                                       true_cell_cluster_allocation, #<- generated_data$true_cell_clust
                                       # true_target_gene_allocation, #  <- generated_data$true_target_gene_allocation
                                       # correct_clustering=correct_clustering,
                                       seed=seeds[i_seed],
                                       do_biclust_with_regulators = do_biclust_with_regulators,
                                       include_regulators_in_results = include_regulators_in_results),
      error = function(e) {
        warning(paste0("Error in bcplaid for c_seed=", i_seed, e$message))
        return(NULL)
      }
    )
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



#
#
# # Find acceptable hyperparams -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# find_hyperparams <- FALSE
# if(find_hyperparams){
#   row.release_vals      = seq(from  = 0, to = 1, length.out = 10)
#   col.release_vals      = seq(from  = 0, to = 1, length.out = 10)
#
#   tuning_ris <- matrix(0,length(row.release_vals), length(col.release_vals))
#
#   iter = 0
#
#   for(rowparam in seq_along(row.release_vals)){
#     for(colparam in seq_along(col.release_vals)){
#
#       iter = iter + 1
#
#       print(paste0("RUNNING BCPLAID ITERATION ", iter, " OF ", length(row.release_vals) * length(col.release_vals), "."))
#
#       biclustbiclust <- function(data){ #centralise this function call so that we only need to set arguments once
#         biclust::biclust(
#           data,
#           method           = BCPlaid(),
#           cluster          ="b",
#           fit.model        = y ~ m + a + b,
#           background       = TRUE,
#           background.layer = NA,
#           background.df    = 1,
#           # row.release      = 0.1,
#           # col.release      = 0.2,
#           row.release      = row.release_vals[rowparam],
#           col.release      = col.release_vals[colparam],
#           shuffle          = 2,
#           back.fit         = 2,
#           max.layers       = 50,
#           iter.startup     = 5,
#           iter.layer       = 100,
#           verbose          = FALSE
#         )
#       }
#       i_seed <- 1234
#
#       tuning_ris  <-tryCatch(
#         expr =  get_stats_biclustbiclust(
#           biclust_input_data     = t(d),
#           n_target_genes     = length(all_ares[[first_converged]][[1]]$call$ind_targetgenes),
#           ind_targetgenes    = all_res[[first_converged]][[1]]$call$ind_targetgenes,
#           n_total_cells      = length(all_res[[first_converged]][[1]]$call$cell_id),
#           n_target_gene_clusters = all_res[[first_converged]][[1]]$call$n_target_gene_clusters,
#           n_cell_clusters        = all_res[[first_converged]][[1]]$call$n_cell_clusters,
#           true_cell_cluster_allocation = all_res[[first_converged]][[1]]$call$true_cell_cluster_allocation,
#           # true_cell_cluster_allocation, #<- generated_data$true_cell_clust
#           # true_target_gene_allocation, #  <- generated_data$true_target_gene_allocation
#           # correct_clustering     = scenarios[[10]]$correct_clustering,
#           seed=i_seed,
#           do_biclust_with_regulators = TRUE,
#           include_regulators_in_results = FALSE
#         )$RI_cell_clustering_biclustbiclust,
#         error = function(e) {
#           warning(paste0("Error in bcplaid for c_seed=", i_seed, e$message))
#           return(NULL)
#         }
#       )
#     }
#   }
# }

#select the best parameter pair
#
# biclustbiclust <- function(data){ #centralise this function call so that we only need to set arguments once
#   biclust::biclust(
#     data,
#     method           = BCPlaid(),
#     cluster          ="b",
#     fit.model        = y ~ m + a + b,
#     background       = TRUE,
#     background.layer = NA,
#     background.df    = 1,
#     # row.release      = 0.1,
#     # col.release      = 0.2,
#     row.release      = 0.3,
#     col.release      = 0.3,
#     shuffle          = 2,
#     back.fit         = 2,
#     max.layers       = 50,
#     iter.startup     = 5,
#     iter.layer       = 100,
#     verbose          = FALSE
#   )
# }



# run

raw_printoutput_path_bcplaid <- file.path(local_data, "output_bcplaid_pbmc_300.txt")
all_res_path_bcplaid <- file.path(local_data, "run300_bcplaid.rds")

if(!file.exists(raw_printoutput_path_bcplaid) || !file.exists(all_res_path_bcplaid) || redo_flag){

  sink(raw_printoutput_path_bcplaid, split=TRUE)

  first_converged <- which(converged)[1]

  # # all_res_bcplaid <- get_stats_biclustbiclust_seeds(
  #   biclust_input_data     = d,
  #   n_target_genes     = length(all_res[[first_converged]][[1]]$call$ind_targetgenes),
  #   ind_targetgenes    = all_res[[first_converged]][[1]]$call$ind_targetgenes,
  #   n_total_cells      = length(all_res[[first_converged]][[1]]$call$cell_id),
  #   n_target_gene_clusters = all_res[[first_converged]][[1]]$call$n_target_gene_clusters,
  #   n_cell_clusters        = all_res[[first_converged]][[1]]$call$n_cell_clusters,
  #   true_cell_cluster_allocation = all_res[[first_converged]][[1]]$call$true_cell_cluster_allocation,
  #   # true_cell_cluster_allocation, #<- generated_data$true_cell_clust
  #   # true_target_gene_allocation, #  <- generated_data$true_target_gene_allocation
  #   # correct_clustering     = scenarios[[10]]$correct_clustering,
  #   seeds                   =  seeds,
  #   do_biclust_with_regulators = TRUE,
  #   include_regulators_in_results = FALSE
  # # )
  #

  all_res_bcplaid <- vector(mode = "list", length = length(seeds))

  for(i_seed in seq_along(seeds)){

    print(paste0("RUNNING BCPLAID ITERATION ", i_seed, " OF ", length(seeds), "."))

    all_res_bcplaid[[i_seed]] <-tryCatch(
      expr =  get_stats_biclustbiclust(
        biclust_input_data     = t(d),
        n_target_genes     = length(all_res[[first_converged]][[1]]$call$ind_targetgenes),
        ind_targetgenes    = all_res[[first_converged]][[1]]$call$ind_targetgenes,
        n_total_cells      = length(all_res[[first_converged]][[1]]$call$cell_id),
        n_target_gene_clusters = all_res[[first_converged]][[1]]$call$n_target_gene_clusters,
        n_cell_clusters        = all_res[[first_converged]][[1]]$call$n_cell_clusters,
        true_cell_cluster_allocation = all_res[[first_converged]][[1]]$call$true_cell_cluster_allocation,
        # true_cell_cluster_allocation, #<- generated_data$true_cell_clust
        # true_target_gene_allocation, #  <- generated_data$true_target_gene_allocation
        # correct_clustering     = scenarios[[10]]$correct_clustering,
        seed=seeds[i_seed],
        do_biclust_with_regulators = TRUE,
        include_regulators_in_results = FALSE
      ),
      error = function(e) {
        warning(paste0("Error in bcplaid for c_seed=", i_seed, e$message))
        return(NULL)
      }
    )
  }


  saveRDS(all_res_bcplaid, all_res_path_bcplaid)
  while (sink.number() > 0) {
    sink()
    sink(file = NULL)
  }
}else{
  all_res_bcplaid <- readRDS(all_res_path_bcplaid)
}


RI_values_plaid <- unlist(sapply(all_res_bcplaid, function(x) x$RI_cell_clustering_biclustbiclust))


### make plots




# Plot info -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
svg(file.path(output_path, paste0("pbmc10k_RI_bcplaid.svg")), width=4, height=3)
boxplot(RI_values_plaid ,
        ylab = "RIs",
        xlab = "BCPlaid",
        # main = "Distribution of rand indexes (RIs)",
        col = c("lightblue"))
dev.off()


# plot_heatmap ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
best_results <- all_res_bcplaid[[which(RI_values_plaid==max(RI_values_plaid, na.rm=TRUE))[1]]]

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


constructed_plots <- plot_biclust_heatmap(biclust_results_matrix                = best_results$biclust_results_matrix,
                                          RI_cell_clustering_biclustbiclust     = best_results$RI_cell_clustering_biclustbiclust,
                                          RI_gene_clustering_biclustbiclust_all = best_results$RI_gene_clustering_biclustbiclust_all,
                                          RI_biclust_biclustbiclust             = best_results$RI_biclust_biclustbiclust)

print(constructed_plots)


png(file.path(output_path, paste0("pbmc_BCPlaid_heatmap.png")), width = 100*aspect_ratio, height = 200, units = "px")
print(constructed_plot)
dev.off()
pdf(file.path(output_path, paste0("pbmc_BCPlaid_heatmap.pdf")), width =  0.95*aspect_ratio/2, height = 1.9)
print(constructed_plot)
dev.off()
