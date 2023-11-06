#!/usr/bin/Rscript

#' Randomise cluster labels
#'
#' Take a vector of integer cluster labels from 1 to x, not skipping any numbers.
#' Assign a fraction fraction_randomised of them a different cluster label.
#' @param cluster_labels Vector of cluster labels, integers, 1 to x, not skipping any numbers.
#' @param fraction_randomised The fraction of each cluster you want reassigned.
#' @return A vector with randomized cluster labels
#' @examples
#' cluster_labels <- rep(c(1,2,3,4), 10)
#' res <- randomise_cluster_labels(cluster_labels = cluster_labels,
#'                                 fraction_randomised  = 0.5)
#' print(cluster_labels)
#' print(res)
#' @export
randomise_cluster_labels <- function(cluster_labels,  # dat$true_cell_cluster_allocation
                                     fraction_randomised = 0.05) {

  disturbed_initial_cell_clust <- cluster_labels
  n_cell_clusters <- length(unique(cluster_labels))

  for (i_cluster in 1:n_cell_clusters) {
    indexes_of_cluster <- which(cluster_labels == i_cluster)
    some_of_those_indexes <- sample(indexes_of_cluster,
                                    size = as.integer(length(indexes_of_cluster) * fraction_randomised),
                                    replace = F)
    disturbed_initial_cell_clust[some_of_those_indexes] <- sample((1:n_cell_clusters)[-i_cluster],  # This makes a vector with all cell cluster labels removing the current one
                                                                  size = length(some_of_those_indexes),
                                                                  replace = T)
  }
  return(disturbed_initial_cell_clust)
}


# runs only when script is run by itself
# || interactive()
if (sys.nframe() == 0){
  # ... do main stuff
  cluster_labels <- rep(c(1,2,3,4), 10)
  res <- randomise_cluster_labels(cluster_labels = cluster_labels,
                                  fraction_randomised  = 0.5)
  print(cluster_labels)
  print(res)
}