R_path <- here::here("R")
source(paste0(R_path,"/generate_dummy_data_for_scregclust.R"))
library(plyr)

generate_dummy_data_for_cell_clustering <- function(
    n_cell_clusters = 3,
    n_target_gene_clusters = c(3,4,5),  # Number of target gene clusters in each cell cluster
    n_target_genes = 40,
    n_regulator_genes = 20,
    n_cells = c(1000,5000,10000),
    regulator_means = c(1,2,3),  # For generating dummy data, regulator mean in each cell cluster
    coefficient_means = list(c(1,20,30), c(1,2,3,4), c(1,2,3,4,5)),  # For generating dummy data, coefficient means in each cell cluster
    disturbed_fraction=0.1  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
){

  # Set variables ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  true_cluster_allocation = rep(1:n_cell_clusters, times=n_cells)
  total_n_cells = sum(n_cells)

  # Generate dummy data for each cell cluster that we want ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  dummy_data <- vector(mode = "list", length = n_cell_clusters)
  betas  <- dummy_data
  for(i_cluster in 1:n_cell_clusters){
    print(paste("Generating data for cell cluster", i_cluster))
    dummy_data[[i_cluster]] <- generate_dummy_data_for_scregclust(n_target_genes,
                                                                  n_regulator_genes,
                                                                  n_cells                = n_cells[i_cluster],
                                                                  n_target_gene_clusters = n_target_gene_clusters[i_cluster],
                                                                  regulator_mean         = regulator_means[i_cluster],
                                                                  coefficient_mean       = coefficient_means[[i_cluster]])
    # list of list of models. with dimension regulators x targets
    # Top level list is each cell cluster
      # in each of those elements is the model generating each target gene cluster in that cell cluster
    betas[[i_cluster]]  <- dummy_data[[i_cluster]]$B
  }

  names(dummy_data[[1]])

  dim(dummy_data[[1]]$B)

  # Create Z_r and Z_t from dummy data
  Z_t <- dummy_data[[1]]$Z_t
  Z_r <- dummy_data[[1]]$Z_r

  #
  scregclust::scregclust(
    expression = rbind(t(Z_t), t(Z_r)),    #scRegClust wants this form
    genesymbols = 1:(n_target_genes+n_regulator_genes),               #gene row numbers
    is_regulator = (1:(n_target_genes+n_regulator_genes) > n_target_genes) + 0, #vector indicating which genes are regulators
    n_cl        = n_target_gene_clusters,
    penalization = 0.001,
    verbose = TRUE
  )-> scRegOut

  str(scRegOut)
  str(scRegOut$results)
  str(scRegOut$results[[1]])
  str(scRegOut$results[[1]]$output)
  scRegOut$results[[1]]$output[[1]]$coeffs
  scRegOut$results[[1]]$output[[1]]$models
  scRegOut$results[[1]]$output[[1]]$cluster[!is_regulator]
  betas_screg <-

  if (n_cell_clusters>1){
    for(i_cluster in 2:n_cell_clusters){
      Z_t <- rbind(Z_t, dummy_data[[i_cluster]]$Z_t)
      Z_r <- rbind(Z_r, dummy_data[[i_cluster]]$Z_r)
    }
  }
  dat <- cbind(Z_t, Z_r)

  # Split into train and test data for cell clustering
  # Skip for now
  # cell_data_split    <- sample(c(1,2), nrow(Z_t), replace = T)
  # train_indices      <- which(cell_data_split == 1)
  # train_dat          <- dat[,train_indices]

  # Get initial cell clustering
  initial_cell_clust <- kmeans(dat, n_cell_clusters)$cluster

  # initial_cell_clust <- sample(1:n_cell_clusters, n_cells, replace = T)


  # Kod för att flytta 1% av cellerna i varje kluster till ett annat kluster.
  disturbed_initial_cell_clust <- initial_cell_clust
  if(disturbed_fraction>0){
    for(i_cluster in 1:n_cell_clusters){
      indexes_of_cluster <- which(initial_cell_clust == i_cluster)
      some_of_those_indexes <- sample(indexes_of_cluster, size=as.integer(length(indexes_of_cluster)*disturbed_fraction), replace = F)
      disturbed_initial_cell_clust[some_of_those_indexes] <- sample(c(1:n_cell_clusters)[-i_cluster], size=length(some_of_those_indexes), replace=T)
    }
  }

  return(list(disturbed_initial_cell_clust=disturbed_initial_cell_clust,
              initial_cell_clust=initial_cell_clust,
              true_cell_clust = true_cell_clust,
              true_betas =  betas,
              dat=dat))
}
