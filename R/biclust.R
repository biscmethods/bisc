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

#' Biclust
#'
#' Cluster cells based with Classification EM based on lm with penality
#' @param dat
#' @param max_iter
#' @param initial_clustering
#' @param n_target_genes
#' @param n_regulator_genes
#' @param n_total_cells
#' @param n_cell_clusters
#' @param ind_targetgenes
#' @param ind_reggenes
#' @param output_path
#' @param penalization_lambda
#' @return Nothing yet, maybe cluster labels
#' @examples
#' NULL
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
  cell_cluster_history[, 'Cell ID'] <- seq_along(initial_clustering)  # Set cell names
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

    for (cell_cluster in 1:n_cell_clusters) {
      cell_cluster_rows <- which(current_cell_cluster_allocation == cell_cluster)
      cell_cluster_target_genes <- as.matrix(dat[cell_cluster_rows, ind_targetgenes])
      cell_cluster_regulator_genes <- as.matrix(dat[cell_cluster_rows, ind_reggenes])
      models[[cell_cluster]] <- lm(formula = 'cell_cluster_target_genes ~ 0 + cell_cluster_regulator_genes',
                                   data = environment())
    }

    # For all cells, calculate the likelihood of coming from the model corresponding to each
    likelihood <- matrix(data = 0, nrow = nrow(dat), ncol = n_cell_clusters)


    # Calculate the residual target gene variance for each gene and cluster
    # (so just one gene).

    # Pre-allocate residual variance estimates
    target_genes_residual_var <- matrix(data = 0, nrow = n_cell_clusters, ncol = n_target_genes)
    # dat is dat <- cbind(target_expression, regulator_expression), e.g. a 2x100, with e.g. the first 50 rows being true cell cluster 1
    # 100x2 * 2x1

    # TODO: This calculates one variance value for each and every target gene type for every cell cluster. Is that correct?
    # Output: n_cell_clusters x n_target_genes
    for (cell_cluster in 1:n_cell_clusters) {
      current_rows <- which(current_cell_cluster_allocation == cell_cluster)
      current_regulator_genes <- as.matrix(dat[current_rows, ind_reggenes])
      current_target_genes <- as.matrix(dat[current_rows, ind_targetgenes])
      cell_cluster_betas <- models[[cell_cluster]]$coefficients

      predicted_values <- current_regulator_genes %*% cell_cluster_betas

      residuals <- current_target_genes - predicted_values

      target_genes_residual_var[cell_cluster,] <- diag(var(residuals))
    }

    # Now to actually calculate predicted or 'predicted' r2
    for (cell in seq_len(nrow(dat))) {
      for (cell_cluster in seq_len(n_cell_clusters)) {
        cell_regulator_genes <- as.matrix(dat[cell, ind_reggenes])
        cell_cluster_betas <- models[[cell_cluster]]$coefficients
        observed_value <- as.matrix(dat[cell, ind_targetgenes])
        cell_cluster_target_genes_residual_var <- target_genes_residual_var[cell_cluster, , drop = FALSE]
        # Bug fix hack: remove NA coefficients
        # if (any(is.na(current_betas))) {
        #   NA_coeffs <- unname(which(is.na(current_betas)))
        #   S_ERR <- (dat[cell, ind_targetgenes] - as.vector(c(1, dat[cell, c(-1, -NA_coeffs)])) %*% current_betas[-NA_coeffs])^2
        # }

        predicted_value <- cell_regulator_genes %*% cell_cluster_betas

        cell_cluster_betas_vector_1norm <- sum(abs(cell_cluster_betas))

        cell_squared_error <- (observed_value - predicted_value)^2

        penalization <- penalization_lambda * cell_cluster_betas_vector_1norm / cell_cluster_target_genes_residual_var

        # TODO: Figure out what the formula should be, and make sure you have which.min or which.max that is correct later around line 220.
        # likelihood[cell,cell_cluster] <- squared_error / 2 / target_genes_residual_var[cell_cluster] - penalization #negative penalty as higher likelihood is better

        # Here we are optimizing the penalized NEGATIVE likelyhood, so penalty is positive
        temp_likelihood <- as.numeric(log(cell_cluster_target_genes_residual_var) / 2 +
                                        cell_squared_error / (2 * target_genes_residual_var[cell_cluster]) +
                                        penalization)
        likelihood[cell, cell_cluster] <- sum(temp_likelihood)

      }
    }

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

    print(paste("Iteration", i_main))
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


  cell_cluster_history_plotting <- cbind('Cell ID' = cell_cluster_history[, 1],
                                         dat$true_cell_cluster_allocation,
                                         cell_cluster_history[, c(2, 3, 4)])
  png(file.path(output_path, paste0("Alluvial_diagram_lambda_", round(penalization_lambda, 3), ".png")))
  plot_cluster_history(cell_cluster_history = cell_cluster_history_plotting)
  dev.off()
}
