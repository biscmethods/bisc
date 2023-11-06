#!/usr/bin/Rscript

library(tidyverse)
library(aricode)  # To calculate rand index
library(ggplot2)  # To plot things #TODO: What things?
library(ggalluvial)  # To plot thingsv #TODO: What things?
library(reshape2)
library(here)  # To work with paths
library(ggfortify)  # For pca-plot

# Get absolute path where script is located, by using relative paths.
R_path <- here::here("R")
source(file.path(R_path, "plot_cluster_history.R"))
source(file.path(R_path, "plot_loglikelihood.R"))


# Compare two matrixes of floats
matequal <- function(x, y, tol = 1e-5) {
  both_matrices_same_dim <- is.matrix(x) &&
    is.matrix(y) &&
    all(dim(x) == dim(y))
  if (both_matrices_same_dim) {
    m <- abs(x - y) <= tol
    all_equal_floats <- all(m)
    print(paste("Number of equal floats:", sum(m), "/", length(m)))
    print("Start of first matrix")
    print(x[1:10, 1:3])
    print("Start of second matrix")
    print(y[1:10, 1:3])
    if (all_equal_floats) {
      print("Matrices are equal.")
      return(TRUE)
    }
  }else {
    print("Matrices are not equal.")
  }
  return(FALSE)
}

loglikelihood_calc <- function(dat,
                               models,
                               target_genes_residual_var,
                               penalization_lambda,
                               n_cell_clusters,
                               ind_reggenes,
                               ind_targetgenes) {
  # For all cells, calculate the likelihood of coming from the model corresponding to each
  likelihood <- matrix(data = 0, nrow = nrow(dat), ncol = n_cell_clusters)

  for (i_cell in seq_len(nrow(dat))) {
    for (i_cell_cluster in seq_len(n_cell_clusters)) {
      cell_regulator_genes <- as.matrix(dat[i_cell, ind_reggenes])  # 1xr
      cell_cluster_betas <- models[[i_cell_cluster]]$coefficients  # rxt
      observed_value <- as.matrix(dat[i_cell, ind_targetgenes]) # 1xt
      cell_cluster_target_genes_residual_var <- target_genes_residual_var[i_cell_cluster, , drop = FALSE]  # 1xt
      predicted_value <- cell_regulator_genes %*% cell_cluster_betas  # 1xt
      cell_cluster_betas_vector_1norm <- sum(abs(cell_cluster_betas))  # 1x1
      cell_squared_error <- (observed_value - predicted_value)^2  # 1xt
      penalization <- penalization_lambda * cell_cluster_betas_vector_1norm / cell_cluster_target_genes_residual_var  # 1xt
      # Here we are optimizing the penalized NEGATIVE likelyhood, so penalty is positive

      temp_likelihood <- as.numeric(log(cell_cluster_target_genes_residual_var) / 2 +
                                      cell_squared_error / (2 * cell_cluster_target_genes_residual_var) +
                                      penalization)  # 1xt
      # print("")
      # print("")
      # print("")
      # print("")
      # print("")
      # print("")
      # print("cell_regulator_genes")
      # print(str(cell_regulator_genes))
      # print("cell_cluster_betas")
      # print(str(cell_cluster_betas))
      # print("observed_value")
      # print(str(observed_value))
      # print("cell_cluster_target_genes_residual_var")
      # print(str(cell_cluster_target_genes_residual_var))
      # print("predicted_value")
      # print(str(predicted_value))
      # print("cell_cluster_betas_vector_1norm")
      # print(str(cell_cluster_betas_vector_1norm))
      # print("cell_squared_error")
      # print(cell_squared_error)
      # print("penalization")
      # print(str(penalization))
      # print("temp_likelihood")
      # print(str(temp_likelihood))
      # stop("STOP")
      likelihood[i_cell, i_cell_cluster] <- sum(temp_likelihood)

    }
  }
  return(likelihood)
}

loglikelihood_calc_matrix <- function(dat,
                                      models,
                                      target_genes_residual_var,
                                      penalization_lambda,
                                      n_cell_clusters,
                                      ind_reggenes,
                                      ind_targetgenes) {
  # For all cells, calculate the likelihood of coming from the model corresponding to each
  n_cells <- nrow(dat)
  likelihood <- matrix(data = 0, nrow = n_cells, ncol = n_cell_clusters)
  for (i_cell_cluster in seq_len(n_cell_clusters)) {
    cell_regulator_genes <- as.matrix(dat[, ind_reggenes])  # 1x20 -> c x 20
    cell_cluster_betas <- models[[i_cell_cluster]]$coefficients  # 20x4
    observed_value <- as.matrix(dat[, ind_targetgenes])  # 1x4 -> cx4
    cell_cluster_target_genes_residual_var <- target_genes_residual_var[i_cell_cluster, , drop = FALSE]  # 1x4
    predicted_value <- cell_regulator_genes %*% cell_cluster_betas  # 1x4 -> cx4
    cell_cluster_betas_vector_1norm <- sum(abs(cell_cluster_betas))  # 1x1
    cell_squared_error <- (observed_value - predicted_value)^2  # 1x4 -> cx4
    penalization <- penalization_lambda * cell_cluster_betas_vector_1norm / cell_cluster_target_genes_residual_var  # 1x4
    # Extrude matrix so we can divide them by element and such
    cell_cluster_target_genes_residual_var <- do.call(rbind, replicate(n_cells, cell_cluster_target_genes_residual_var, simplify = FALSE))  # cx4
    penalization <- do.call(rbind, replicate(n_cells, penalization, simplify = FALSE))  #  cx4

    # Here we are optimizing the penalized NEGATIVE likelyhood, so penalty is positive
    temp_likelihood <- log(cell_cluster_target_genes_residual_var) / 2 +
      cell_squared_error / (2 * cell_cluster_target_genes_residual_var) +
      penalization  # 1x4 -> cx4
    # print("")
    # print("")
    # print("")
    # print("")
    # print("")
    # print("")
    # print("cell_regulator_genes")
    # print(str(cell_regulator_genes))
    # print("cell_cluster_betas")
    # print(str(cell_cluster_betas))
    # print("observed_value")
    # print(str(observed_value))
    # print("cell_cluster_target_genes_residual_var")
    # print(str(cell_cluster_target_genes_residual_var))
    # print("predicted_value")
    # print(str(predicted_value))
    # print("cell_cluster_betas_vector_1norm")
    # print(str(cell_cluster_betas_vector_1norm))
    # print("cell_squared_error")
    # print(str(cell_squared_error))
    # print("penalization")
    # print(str(penalization))
    # print("target_genes_residual_var")
    # print(str(target_genes_residual_var))
    # print("log(cell_cluster_target_genes_residual_var)/2")
    # print(log(cell_cluster_target_genes_residual_var)/2)
    # print("temp_likelihood")
    # print(str(temp_likelihood))
    # stop("STOP")

    likelihood[, i_cell_cluster] <- rowSums(temp_likelihood)

  }

  return(likelihood)
}

#' Biclust
#'
#' Cluster cells based with Classification EM based on lm with penality
#' @param dat Main data
#' @param max_iter Max number of iterations to run
#' @param initial_clustering The start clustering
#' @param n_target_genes Number of target genes
#' @param n_regulator_genes Number of regulator genes
#' @param n_total_cells Number of total cells
#' @param n_cell_clusters Number of cell clusters
#' @param ind_targetgenes Indexes of where the target gene columns in dat are
#' @param ind_reggenes Indexes of where the regulator gene columns in dat are
#' @param output_path Output path of plots
#' @param penalization_lambda The penalization lambda
#' @return Nothing yet, maybe cluster labels
#' @export
biclust <- function(dat = dat,
                    max_iter = 50,
                    initial_clustering,
                    n_target_genes,
                    n_regulator_genes,
                    n_total_cells,
                    n_cell_clusters,
                    ind_targetgenes,
                    ind_reggenes,
                    output_path,
                    penalization_lambda = 0.5) {

  # Preallocate cluster history
  cell_cluster_history <- tibble::tibble(dat$cell_id, dat$true_cell_cluster_allocation, initial_clustering)
  colnames(cell_cluster_history) <- c("cell_id", "True allocation", "Disturbed allocation")
  initial_column_padding <- ncol(as.matrix(initial_clustering)) + 1  # +1 Because we have an index column that is not an index column it's an ID column
  cell_cluster_history <- data.frame(matrix(NA, nrow = nrow(as.matrix(initial_clustering)), ncol = max_iter + initial_column_padding))
  colnames(cell_cluster_history) <- c("Cell ID", 'Initial clustering', paste0("Iteration ", 1:max_iter))
  cell_cluster_history[, 'Cell ID'] <- seq_along(initial_clustering) # Set cell names
  cell_cluster_history[, 'Initial clustering'] <- initial_clustering
  cell_cluster_history <- tibble::as_tibble(cell_cluster_history)


  # Pre-allocate all r2 matrices for later analysis if feasible
  likelihood_all <- vector(mode = "list", length = max_iter)

  # Set the current cell clustering
  current_cell_cluster_allocation <- initial_clustering

  stop_iterating_flag <- 0  # Flag if we have converged
  for (i_main in 1:max_iter) {

    # Fit model to each cell cluster
    models <- vector(mode = "list", length = n_cell_clusters)

    for (i_cell_cluster in 1:n_cell_clusters) {
      cell_cluster_rows <- which(current_cell_cluster_allocation == i_cell_cluster)
      cell_cluster_target_genes <- as.matrix(dat[cell_cluster_rows, ind_targetgenes])
      cell_cluster_regulator_genes <- as.matrix(dat[cell_cluster_rows, ind_reggenes])
      models[[i_cell_cluster]] <- lm(formula = 'cell_cluster_target_genes ~ 0 + cell_cluster_regulator_genes',
                                     data = environment())
    }


    # Calculate the residual target gene variance for each gene and cluster
    # (so just one gene).

    # Pre-allocate residual variance estimates
    target_genes_residual_var <- matrix(data = 0, nrow = n_cell_clusters, ncol = n_target_genes)
    # dat is dat <- cbind(target_expression, regulator_expression), e.g. a 2x100, with e.g. the first 50 rows being true cell cluster 1
    # 100x2 * 2x1

    # TODO: This calculates one variance value for each and every target gene type for every cell cluster. Is that correct?
    # Output: n_cell_clusters x n_target_genes
    start.time <- Sys.time()
    for (i_cell_cluster in seq_len(n_cell_clusters)) {
      current_rows <- which(current_cell_cluster_allocation == i_cell_cluster)
      current_regulator_genes <- as.matrix(dat[current_rows, ind_reggenes])
      current_target_genes <- as.matrix(dat[current_rows, ind_targetgenes])
      cell_cluster_betas <- models[[i_cell_cluster]]$coefficients

      predicted_values <- current_regulator_genes %*% cell_cluster_betas

      residuals <- current_target_genes - predicted_values

      target_genes_residual_var[i_cell_cluster,] <- diag(var(residuals))
    }
    time_taken <- round(Sys.time() - start.time, 2)
    print(paste("Iteration", i_main, "res var", time_taken))

    # Now to actually calculate predicted or 'predicted' r2
    # start.time <- Sys.time()
    # likelihood <- loglikelihood_calc(dat,
    #                                  models,
    #                                  target_genes_residual_var,
    #                                  penalization_lambda,
    #                                  n_cell_clusters,
    #                                  ind_reggenes,
    #                                  ind_targetgenes)
    # time_taken <- round(Sys.time() - start.time, 2)
    # print(paste("Iteration", i_main, "likelihood", time_taken))

    start.time <- Sys.time()
    likelihood <- loglikelihood_calc_matrix(dat,
                                                   models,
                                                   target_genes_residual_var,
                                                   penalization_lambda,
                                                   n_cell_clusters,
                                                   ind_reggenes,
                                                   ind_targetgenes)
    time_taken <- round(Sys.time() - start.time, 2)
    print(paste("Iteration", i_main, "likelihood matrix form", time_taken))

    # matequal(likelihood, likelihood_matrix)

    likelihood_all[[i_main]] <- likelihood

    # r2plot(iteration = i_main,
    #      r2 = likelihood,
    #      prev_cell_clust = dat$true_cell_cluster_allocation)

    # scatter_plot_loglikelihood()
    # hist_plot_loglikelihood()

    # Update cluster allocations
    updated_cell_clust <- sapply(seq_len(nrow(likelihood)), function(row) which.min(likelihood[row,]))

    # Update data in cell_cluster_history
    cell_cluster_history[, i_main + initial_column_padding] <- updated_cell_clust
    aricode::RI(dat$true_cell_cluster_allocation, updated_cell_clust)

    # Check convergence of cluster labels
    # Compare clusters with with previous iterations so we can exit if we seen this allocation before
    for (prev_clustering in ((i_main - 1):0)) {
      print(paste0('Comparing with iteration ', prev_clustering))
      rand_index <- aricode::RI(updated_cell_clust,
                                as.matrix(cell_cluster_history)[, prev_clustering + initial_column_padding]
      )
      if (rand_index == 1) {
        print("Cell clustering from iteration same as some previous iteration. Exiting.")
        print(paste0("RI of ", rand_index,
                     " when comparing iteration ", i_main,
                     " to iteration ", prev_clustering))
        stop_iterating_flag <- T
        break
      }
    }

    if (stop_iterating_flag) {
      # Clean up cluster history
      cell_cluster_history <- cell_cluster_history[, colSums(is.na(cell_cluster_history)) == 0, drop = FALSE]
      # Stop iterations/exit function
      break
    }

    scatter_plot_loglikelihood(dat,
                               likelihood,
                               n_cell_clusters,
                               penalization_lambda,
                               output_path,
                               i_main)
    hist_plot_loglikelihood(dat,
                            likelihood,
                            n_cell_clusters,
                            penalization_lambda,
                            output_path,
                            i_main)
  }

  start.time <- Sys.time()
  cell_cluster_history_plotting <- cbind('Cell ID' = cell_cluster_history[, 1],
                                         dat$true_cell_cluster_allocation,
                                         cell_cluster_history[, c(2, 3, 4)])
  png(file.path(output_path, paste0("Alluvial_diagram_lambda_", round(penalization_lambda, 3), ".png")))
  plot_cluster_history(cell_cluster_history = cell_cluster_history_plotting)
  dev.off()
  time_taken <- round(Sys.time() - start.time, 2)
  print(paste("Iterations complete, Alluvial plot", time_taken))
}
