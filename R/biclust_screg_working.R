#!/usr/bin/Rscript

library(tidyverse)
library(aricode)     # To calculate rand index
library(ggplot2)     # To plot things #TODO: What things? literally anything
library(ggalluvial)  # To plot thingsv #TODO: What things?
library(reshape2)
library(here)        # To work with paths
library(ggfortify)   # For pca-plot
library(pracma)      # For pseudo inverse
library(stats)       # For kmeans
library(ppclust)
library(clusterSim)  # for db

# Get absolute path where script is located, by using relative paths.
R_path <- here::here("R")
source(file.path(R_path, "plot_cluster_history.R"))
source(file.path(R_path, "plot_loglikelihood.R"))
source(file.path(R_path, "generate_dummy_data_for_cell_clustering.R"))


# Compare two matrixes of floats
matequal <- function(x, y, tol = 1e-5) {
  both_matrices_same_dim <- is.matrix(x) &&
    is.matrix(y) &&
    all(dim(x) == dim(y))
  if (both_matrices_same_dim) {
    m <- abs(x - y) <= tol
    all_equal_floats <- all(m)
    print(paste(" Number of equal floats:", sum(m), "/", length(m)), quote = FALSE)
    print(" Start of first matrix", quote = FALSE)
    print(x[1:10, 1:3], quote = FALSE)
    print(" Start of second matrix", quote = FALSE)
    print(y[1:10, 1:3], quote = FALSE)
    if (all_equal_floats) {
      print(" Matrices are equal.", quote = FALSE)
      return(TRUE)
    }
  }else {
    print(" Matrices are not equal.", quote = FALSE)
  }
  return(FALSE)
}

likelihood_calc <- function(dat,
                            models,
                            target_genes_residual_var,
                            penalization_lambda,
                            n_cell_clusters,
                            ind_reggenes,
                            ind_targetgenes) {

  # Actually calculates loglikelihoods at the moment

  # For all cells, calculate the likelihood of coming from the model corresponding to each
  likelihood <- matrix(data = 0, nrow = nrow(dat), ncol = n_cell_clusters)

  for (i_cell in seq_len(nrow(dat))) {
    for (i_cell_cluster in seq_len(n_cell_clusters)) {

      # Get regulators and targets and coefficients
      cell_regulator_genes <- as.matrix(dat[i_cell, ind_reggenes, drop = FALSE])  # 1xr
      observed_value <- as.matrix(dat[i_cell, ind_targetgenes, drop = FALSE]) # 1xt
      cell_cluster_betas <- models[[i_cell_cluster]]   # rxt

      # Get residual variances and predicted values
      cell_cluster_target_genes_residual_var <- target_genes_residual_var[i_cell_cluster, , drop = FALSE]  # 1xt
      predicted_value <- cell_regulator_genes %*% cell_cluster_betas  # 1xt

      # Cell_cluster_betas_vector_1norm <- sum(abs(cell_cluster_betas))  # 1x1
      cell_squared_error <- (observed_value - predicted_value)^2  # 1xt
      # penalization <- penalization_lambda * cell_cluster_betas_vector_1norm / cell_cluster_target_genes_residual_var  # 1xt

      # Here we are optimizing the penalized NEGATIVE likelyhood, so penalty is positive
      # temp_likelihood <- as.numeric(log(cell_cluster_target_genes_residual_var) / 2 +
      #                               cell_squared_error / (2 * cell_cluster_target_genes_residual_var) +
      #                               penalization)  # 1xt
      # temp_likelihood <- sum(temp_likelihood)

      # Calcualte probability that each observation comes from each cluster given the current model

      # should triple-check this math

      term1 <- log(2 * pi)                                                         # constant
      term2 <- log(cell_cluster_target_genes_residual_var) / 2                   # vector of length t
      term3 <- cell_squared_error / (cell_cluster_target_genes_residual_var * 2) # vector of length t


      # Sum up and exponentiate back. Add/remove exp for likelihoods/loglikelihoods
      temp_likelihood <- (sum(-term1 - term2 - term3))  # Scalar

      # Actually loglikelihood
      likelihood[i_cell, i_cell_cluster] <- temp_likelihood
    }
  }
  return(likelihood)
}

loglikelihood_calc_matrix <- function(dat,
                                      models,
                                      target_genes_residual_var,
                                      n_cell_clusters,
                                      ind_reggenes,
                                      ind_targetgenes) {
  # For all cells, calculate the likelihood of coming from the model corresponding to each
  n_cells <- nrow(dat)
  loglikelihood <- matrix(data = 0, nrow = n_cells, ncol = n_cell_clusters)
  for (i_cell_cluster in seq_len(n_cell_clusters)) {
    cell_regulator_genes <- as.matrix(dat[, ind_reggenes, drop = FALSE])  # 1xr -> cxr
    observed_value <- as.matrix(dat[, ind_targetgenes, drop = FALSE])  # 1xt -> cxt
    cell_cluster_betas <- models[[i_cell_cluster]] # rxt

    predicted_value <- cell_regulator_genes %*% cell_cluster_betas  # 1xt -> cxt
    cell_squared_error <- (observed_value - predicted_value)^2  # 1xt -> cxt
    cell_squared_error[is.na(cell_squared_error)] <- mean(cell_squared_error, na.rm = TRUE)

    # Extrude 1xt to cxt so we can divide element wise
    cell_cluster_target_genes_residual_var <- target_genes_residual_var[i_cell_cluster, , drop = FALSE]  # 1xt
    cell_cluster_target_genes_residual_var <- do.call(rbind, replicate(n_cells, cell_cluster_target_genes_residual_var, simplify = FALSE))  # cxt
    term1 <- log(2 * pi)                                                         # scalar
    term2 <- log(cell_cluster_target_genes_residual_var) / 2                     # 1xt -> cxt
    term3 <- cell_squared_error / (cell_cluster_target_genes_residual_var * 2)   # 1xt -> cxt

    temp_loglikelihood <- rowSums(-term1 - term2 - term3)  # cx1

    loglikelihood[, i_cell_cluster] <- temp_loglikelihood
  }

  return(loglikelihood)
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
#'
biclust <- function(dat = dat,
                    cell_id,
                    true_cell_cluster_allocation,
                    max_iter = 50,
                    n_target_gene_clusters,
                    initial_clustering,
                    n_target_genes,
                    n_regulator_genes,
                    n_total_cells,
                    n_cell_clusters,
                    ind_targetgenes,
                    ind_reggenes,
                    output_path,
                    penalization_lambda = 0.05,
                    use_weights = TRUE,
                    use_complex_cluster_allocation = FALSE) {

  if (!is.factor(initial_clustering)) {
    stop("The variable initial_clustering needs to be a factor vector.")
  }

  # Preallocate cluster history

  cell_cluster_history <- tibble::tibble(cell_id, true_cell_cluster_allocation, initial_clustering)
  colnames(cell_cluster_history) <- c("cell_id", "True allocation", "Disturbed allocation")
  initial_column_padding <- ncol(as.matrix(initial_clustering)) + 1  # +1 Because we have an index column that is not an index column it's an ID column
  cell_cluster_history <- data.frame(matrix(NA, nrow = nrow(as.matrix(initial_clustering)), ncol = max_iter + initial_column_padding))
  colnames(cell_cluster_history) <- c("Cell ID", 'Initial clustering', paste0("Iteration ", 1:max_iter))
  cell_cluster_history[, 'Cell ID'] <- seq_along(initial_clustering) # Set cell names
  cell_cluster_history[, 'Initial clustering'] <- initial_clustering
  cell_cluster_history <- tibble::as_tibble(cell_cluster_history)


  # Pre-allocate all r2 matrices for later analysis if feasible
  likelihood_all <- vector(mode = "list", length = max_iter)
  weights_all <- vector(mode = "list", length = max_iter)
  BIC_all <- vector(mode = "list", length = max_iter)
  target_genes_residual_var_all <- vector(mode = "list", length = max_iter)
  # Set the current cell clustering
  current_cell_cluster_allocation <- initial_clustering

  stop_iterating_flag <- 0  # Flag if we have converged
  for (i_main in 1:max_iter) {
    print(paste(" Running iteration", i_main), quote = FALSE)
    start.time <- Sys.time()

    ################################################################
    ##### M-step, compute estimates for \pi_k and model parameters #
    ##### \pi_k are cluster proportions                            #
    ##### model parameters are betas and signas from SCREG         #
    ################################################################

    ###### Fit model to each cell cluster ####
    ###### M.1                           ####
    models <- vector(mode = "list", length = n_cell_clusters)

    # #dev
    # models2 <- vector(mode = "list", length = n_cell_clusters)
    # models_eval <- vector(mode = "list", length = n_cell_clusters)

    for (i_cell_cluster in 1:n_cell_clusters) {
      print(paste("  Calculating betas for cell cluster", i_cell_cluster), quote = FALSE)

      cell_cluster_rows <- which(current_cell_cluster_allocation == i_cell_cluster)
      cell_cluster_target_genes <- as.matrix(dat[cell_cluster_rows, ind_targetgenes, drop = FALSE])
      cell_cluster_regulator_genes <- as.matrix(dat[cell_cluster_rows, ind_reggenes, drop = FALSE])

      # For SCREG we can skip weights for now as it involves updating the screg function
      #
      #if (use_weights == FALSE || i_main == 1) {
      #  current_weights <- diag(nrow(cell_cluster_regulator_genes))
      #}
      #else {
      # Current weights is a n_cell x n_cell matrix with the weigths on the diagonal
      #  current_weights <- weights_all[[i_main - 1]][cell_cluster_rows]
      # current_weights <- current_weights / sum(current_weights)
      #  current_weights <- diag(current_weights)
      #}

      # Now we fit the models

      # TODO: Replace this with SCREGCLUST calls
      #models[[i_cell_cluster]] <- ridge_weighted_lm_(
      #  X = cell_cluster_regulator_genes,
      #  Y = cell_cluster_target_genes,
      #  W = current_weights,
      #  lambda = penalization_lambda
      #)
      #  scregclust is function that takes expression = p rows of genes and n columns of cells.
      if (nrow(cell_cluster_target_genes) == 0) {
        stop(paste("Number of cells in cell cluster ", i_cell_cluster, "is 0"))
      }
      screg_out <- scregclust::scregclust(
        expression = rbind(t(cell_cluster_target_genes), t(cell_cluster_regulator_genes)),  # scRegClust wants this form
        genesymbols = 1:(n_target_genes + n_regulator_genes),  # Gene row numbers
        is_regulator = (1:(n_target_genes + n_regulator_genes) > n_target_genes) + 0,  # Vector indicating which genes are regulators
        n_cl = n_target_gene_clusters[[i_cell_cluster]],
        penalization = penalization_lambda,
        noise_threshold = 0.00001,
        verbose = FALSE
      )
      screg_out_betas <- do.call(cbind, screg_out$results[[1]]$output[[1]]$coeffs)  # Merge betas into one matrix
      target_gene_cluster_names <- screg_out$results[[1]]$output[[1]]$cluster[1:n_target_genes]
      n_target_gene_cluster_names <- length(unique(target_gene_cluster_names))

      # if (any(target_gene_cluster_names == -1)) {
      #   stop(paste("Number of target genes put into the rubbish cluster are", sum(target_gene_cluster_names == -1), "out of", length(target_gene_cluster_names)))
      # }
      #
      # if (any(is.na(target_gene_cluster_names))) {
      #   stop("For some reason there are NAs in the target gene cluster allocation vector.")
      # }
      #
      # if (n_target_gene_cluster_names != n_target_gene_clusters[[i_cell_cluster]]) {
      #   stop(paste("Number of found clusters by scregclust was", n_target_gene_cluster_names, "but should be", n_target_gene_clusters[[i_cell_cluster]]))
      # }

      # Create cluster_indexes which maps columns in screg_out_betas to correct places under some assumptions of order:
      # target_gene_cluster_names, e.g. 2 2 1 3 1 3,
      # becomes cluster_indexes: 3 4 1 5 2 6
      # CAN ALSO BE DONE WITH THE FOLLOWING CODE THAT WE HAVE USED BEFORE:
      # indexes_of_not_deleted_target_gene_clusters <- which(sapply(out_list[[i_target_cell_cluster]]$results[[1]]$output[[1]]$coeffs, FUN=function(x) !is.null(x)))
      # target_gene_cluster_index <- out_list[[i_target_cell_cluster]]$results[[1]]$output[[1]]$cluster[1:n_target_genes]  # n_target_genes, e.g. 1 1 1 1 3 3 3 3 3 3 3 3 1 3 1 3 3 3 3 3 1 3 1 1 3 1 3 3 1 3
      # target_gene_cluster_index <- unlist(sapply(indexes_of_not_deleted_target_gene_clusters, FUN=function(x) which(target_gene_cluster_index==x)))
      cluster_indexes <- rep(NA, length(target_gene_cluster_names))
      tmp_max_ind <- 0
      for (i_target_gene_cluster in 1:n_target_gene_clusters[[i_cell_cluster]]) {
        ind_tmp <- which(target_gene_cluster_names == i_target_gene_cluster)
        n_target_genes_in_cluster <- length(ind_tmp)
        cluster_indexes[ind_tmp] <- 1:n_target_genes_in_cluster + tmp_max_ind  # Get names of clusters for each target gene
        tmp_max_ind <- n_target_genes_in_cluster
      }

      screg_out_betas <- screg_out_betas[, cluster_indexes]  # Sort the matrix

      # screg_out_sigmas <- do.call(c, screg_out$results[[1]]$output[[1]]$sigmas)
      # screg_out_sigmas <- screg_out_sigmas[cluster_indexes]

      models[[i_cell_cluster]] <- screg_out_betas

      # TODO: alternatively write some predict function and use extracted variance estimates

    }

    #### calculate one variance for each target gene given each model ####
    ###### M.2                                                        ####

    # Pre-allocate residual variance estimates
    target_genes_residual_var <- matrix(data = 0, nrow = n_cell_clusters, ncol = n_target_genes)

    # Output: n_cell_clusters x n_target_genes
    for (i_cell_cluster in seq_len(n_cell_clusters)) {
      print(paste("  Calculating residual variance for cell cluster", i_cell_cluster), quote = FALSE)

      # TODO: here we will need to take SCREGCLUSTs normalisation into account
      current_rows <- which(current_cell_cluster_allocation == i_cell_cluster)
      current_regulator_genes <- as.matrix(dat[current_rows, ind_reggenes, drop = FALSE])
      current_target_genes <- as.matrix(dat[current_rows, ind_targetgenes, drop = FALSE])
      cell_cluster_betas <- models[[i_cell_cluster]] # $coefficients # with our own lm we only save the coeffs anyway

      predicted_values <- current_regulator_genes %*% cell_cluster_betas

      residuals <- current_target_genes - predicted_values

      # target_genes_residual_var[i_cell_cluster,] <- diag(var(residuals))  # maybe not necessary to calculate entire matrix
      # TODO: this should be possible to just extract from SCREG, do that and compare to make sure it's approx same
      tmp_res_var <- colSums(residuals^2) / (length(current_rows) - 1)
      tmp_res_var[is.na(tmp_res_var)] <- mean(tmp_res_var, na.rm = TRUE)  # Make target genes that were thrown in the rubbish cluster have large variance
      target_genes_residual_var[i_cell_cluster,] <- tmp_res_var

      #dev
      # cell_cluster_betas2 <- models2[[i_cell_cluster]]
      # predicted_values2 <- current_regulator_genes %*% cell_cluster_betas2
      # residuals2 <- current_target_genes - predicted_values2
      # target_genes_residual_var2[i_cell_cluster,] <- diag(var(residuals2))
    }
    target_genes_residual_var_all[[i_main]] <- target_genes_residual_var


    # Calculated loglikelihoods
    # M.3

    # likelihood <- likelihood_calc(dat,
    #                               models,
    #                               target_genes_residual_var,
    #                               penalization_lambda,
    #                               n_cell_clusters,
    #                               ind_reggenes,
    #                               ind_targetgenes)

    likelihood <- loglikelihood_calc_matrix(dat,
                                            models,
                                            target_genes_residual_var,
                                            n_cell_clusters,
                                            ind_reggenes,
                                            ind_targetgenes)

    # matequal(likelihood, loglikelihood_matrix)
    cluster_proportions <- vector(length = n_cell_clusters)
    for (i_cell_cluster in seq_len(n_cell_clusters)) {
      current_cluster_proportion <- sum(current_cell_cluster_allocation == i_cell_cluster) / length(current_cell_cluster_allocation)

      # If some cluster has been completely omitted, give it a nonzero proportion anyway for next iteration
      # Even without weights this is good to include for cluster allocation (though it is a bit arbitrary)
      if (current_cluster_proportion == 0) {
        cluster_proportions[i_cell_cluster] <- 0.0001
        stop_iterating_flag <- T
      }else
        cluster_proportions[i_cell_cluster] <- current_cluster_proportion
    }
    # cluster_proportions <- unname(table(current_cell_cluster_allocation) /
    #                                 length(current_cell_cluster_allocation))

    ####################################################################
    ##### E-step #######################################################
    ####### Compute posterior probabilities ############################
    ####### For each cell to belong to each cluster ####################

    # If everything goes well this should work the same for screg as with our
    # lm examples
    weights <- sweep(exp(likelihood), 2, cluster_proportions, "*")

    big_logLhat_in_BIC <- sum(log(rowSums(weights)))
    k_in_BIC <- (n_target_genes + n_regulator_genes) * n_cell_clusters
    n_in_BIC <- n_total_cells
    BIC <- k_in_BIC * log(n_in_BIC) - 2 * big_logLhat_in_BIC

    weights <- sweep(weights, 1, rowSums(weights), "/")

    likelihood <- weights

    weights_all[[i_main]] <- weights
    BIC_all[[i_main]] <- BIC

    likelihood_all[[i_main]] <- likelihood

    if (i_main == 1) {
      no_factor_cluster <- as.numeric(levels(initial_clustering))[initial_clustering]  # This weird line converts a vector with factors to a normal numeric vector.
      # print(str(no_factor_cluster), quote = FALSE)
      # db <- index.DB(likelihood, no_factor_cluster)$DB
      db <- mean(silhouette(no_factor_cluster, dist(likelihood))[, 3])
    }

    ######################################
    ####  Calculate cluster proportions ##
    ######  Step  M.4 ####################

    # cluster_proportions <- unname(table(current_cell_cluster_allocation) /
    #                                 length(current_cell_cluster_allocation))
    #
    # ####################################################################
    # ##### E-step #######################################################
    # ####### Compute posterior probabilities ############################
    # ####### For each cell to belong to each cluster ####################
    #
    # weights <- sweep(exp(likelihood), 2, cluster_proportions, "*")
    # weights <- sweep(weights, 2, colSums(weights), "/")
    # #
    # # # weights_2 <- sweep(likelihood, 2, cluster_proportions, "+")
    # # # weights_2 <- sweep(weights_2, 2, colSums(weights_2), "-")
    # # # print(str(weights_2))
    # # # print(str(weights))
    # #
    # weights_all[[i_main]] <- weights

    ####################################################################
    ##### E-step #######################################################
    ##### update cluster allocations ###################################
    if (use_complex_cluster_allocation) {
      pca_likelihood <- stats::prcomp(likelihood, center = TRUE, scale. = TRUE)$x
      fcm_res <- ppclust::fcm(x = pca_likelihood, centers = n_cell_clusters)
      weights_all[[i_main]] <- fcm_res$u
      updated_cell_clust <- fcm_res$cluster
      # updated_cell_clust <- stats::kmeans(x = pca_likelihood, centers = n_cell_clusters, iter.max = 20, nstart = 50 + n_cell_clusters)$cluster
    }else {
      updated_cell_clust <- sapply(seq_len(nrow(likelihood)),
                                   function(row) which.max(likelihood[row,]))
    }
    updated_cell_clust <- unlist(updated_cell_clust)
    # print("updated_cell_clust")
    # print((str(unlist(updated_cell_clust))))
    #
    # print("current_cell_cluster_allocation")
    # print(str(current_cell_cluster_allocation))

    # Update data in cell_cluster_history

    cell_cluster_history[, i_main + initial_column_padding] <- updated_cell_clust

    # Check convergence of cluster labels
    # Compare clusters with with previous iterations so we can exit if we seen this allocation before
    for (prev_clustering in ((i_main - 1):0)) {
      rand_index <- aricode::RI(updated_cell_clust,
                                as.matrix(cell_cluster_history)[, prev_clustering + initial_column_padding]
      )
      if (rand_index == 1) {
        print(" Cell clustering from iteration same as some previous iteration. Exiting.", quote = FALSE)
        print(paste0(" Rand Index of ", rand_index,
                     " when comparing iteration ", i_main,
                     " to iteration ", prev_clustering), quote = FALSE)
        stop_iterating_flag <- T
        break
      }
    }

    rand_index_true_cluster <- aricode::RI(true_cell_cluster_allocation, updated_cell_clust)
    time_taken <- round(Sys.time() - start.time, 2)
    print(paste0(" Iteration ", i_main, ", took ", time_taken, " seconds", ", Rand index: ", rand_index_true_cluster), quote = FALSE)

    current_cell_cluster_allocation <- as.factor(updated_cell_clust)

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
                               i_main,
                               true_cell_cluster_allocation)
    hist_plot_loglikelihood(dat,
                            likelihood,
                            n_cell_clusters,
                            penalization_lambda,
                            output_path,
                            i_main,
                            true_cell_cluster_allocation)
  }

  start.time <- Sys.time()

  cell_cluster_history_plotting <- cbind('Cell ID' = cell_cluster_history[, 1],
                                         'True cell cluster allocation' = true_cell_cluster_allocation,
                                         cell_cluster_history[, 2:ncol(cell_cluster_history)])
  png(file.path(output_path, paste0("Alluvial_diagram_lambda_", round(penalization_lambda, 3), ".png")),
      width = 1024 + ncol(cell_cluster_history_plotting) * 40, height = 1024, units = "px", res = 150)
  plot_cluster_history(cell_cluster_history = cell_cluster_history_plotting, correct_plot = FALSE)
  dev.off()
  time_taken <- round(Sys.time() - start.time, 2)
  print(paste(" Iterations complete, Alluvial plot took", time_taken, "seconds."), quote = FALSE)

  return(list("rand_index" = rand_index_true_cluster,
              "n_iterations" = i_main,
              "db" = db,
              "BIC" = BIC_all[1:i_main],
              "taget_genes_residual_var" = target_genes_residual_var_all[1:i_main]))
}


# Runs only when script is run by itself
# || interactive()
if (sys.nframe() == 0) {
  # ... do main stuff


  #############################################
  ############ data for dev ###################
  #############################################

  # Set seed for example
  set.seed(3)

  # Set variables ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  n_cell_clusters <- 2
  n_target_gene_clusters <- c(2, 3)  # Number of target gene clusters in each cell cluster
  n_target_genes <- 50
  n_regulator_genes <- 30
  n_cells <- c(10000, 5000)
  regulator_means <- c(5, 1)  # For generating dummy data, regulator mean in each cell cluster
  coefficient_means <- list(c(1, 20), c(5, 20, 100))  # For generating dummy data, coefficient means in each cell cluster
  coefficient_sds <- list(c(0.1, 0.1), c(0.1, 0.1, 0.1))
  true_cluster_allocation <- rep(1:n_cell_clusters, times = n_cells)
  n_total_cells <- sum(n_cells)

  # Generate dummy data for each cell cluster that we want ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  generated_data <- generate_dummy_data_for_cell_clustering(
    n_cell_clusters = n_cell_clusters,
    n_target_gene_clusters = n_target_gene_clusters,  # Number of target gene clusters in each cell cluster
    n_target_genes = n_target_genes,
    n_regulator_genes = n_regulator_genes,
    n_cells = n_cells,
    regulator_means = regulator_means,  # For generating dummy data, regulator mean in each cell cluster
    coefficient_means = coefficient_means,  # For generating dummy data, coefficient means in each cell cluster
    coefficient_sds = coefficient_sds,
    disturbed_fraction = 0.22
  )

  ind_reggenes <- which(c(rep(0, n_target_genes), rep(1, n_regulator_genes)) == 1)
  ind_targetgenes <- which(c(rep(1, n_target_genes), rep(0, n_regulator_genes)) == 1)

  disturbed_initial_cell_clust <- factor(generated_data$disturbed_initial_cell_clust)

  biclust_input_data <- generated_data$dat
  colnames(biclust_input_data) <- c(paste0("r", 1:n_target_genes), paste0("t", 1:n_regulator_genes))
  biclust_input_data <- tibble::as_tibble(biclust_input_data)

  # # These needs to be strings for discrete labels in pca plot
  # data_for_plotting <- tibble::as_tibble(true_cell_cluster_allocation = generated_data$true_cell_clust,
  #                                        biclust_input_data)
  # pca_res <- prcomp(biclust_input_data, scale. = TRUE)
  # p <- ggplot2::autoplot(pca_res, data = data_for_plotting, colour = 'true_cell_cluster_allocation')

  # Set up some variables
  n_cell_clusters <- length(unique(disturbed_initial_cell_clust))
  n_target_genes <- length(ind_targetgenes)
  n_regulator_genes <- length(ind_reggenes)


  #############################################
  ############ end data for dev ###############
  #############################################

  ###########################################
  ####initialise variables for dev ##########
  ##########################################

  max_iter <- 50
  initial_clustering <- disturbed_initial_cell_clust
  n_target_genes <- n_target_genes
  n_regulator_genes <- n_regulator_genes
  n_total_cells <- n_total_cells
  n_cell_clusters <- n_cell_clusters
  ind_targetgenes <- ind_targetgenes
  ind_reggenes <- ind_reggenes
  # output_path <- modded_output_path
  penalization_lambda <- 0.8
  i_cell_cluster <- 1
  i_main <- 1

  use_weights <- TRUE
  use_complex_cluster_allocation <- FALSE

  demo_path <- here::here("demo")
  output_path <- demo_path
  # initial_clustering[1:14900] <- rep(1, 14900)
  # initial_clustering[14901:15000] <- rep(2, 100)
  # print(length(initial_clustering))
  # stop("hej")


  ###########################################
  ###END initialise variables for dev #######
  ###########################################
  biclust_result <- biclust(dat = biclust_input_data,
                            cell_id = 1:n_total_cells,
                            true_cell_cluster_allocation = factor(generated_data$true_cell_clust),
                            max_iter = 50,
                            n_target_gene_clusters,
                            initial_clustering,
                            n_target_genes,
                            n_regulator_genes,
                            n_total_cells,
                            n_cell_clusters,
                            ind_targetgenes,
                            ind_reggenes,
                            output_path,
                            penalization_lambda = penalization_lambda)
  print(paste("rand_index for result vs true cluster:", biclust_result$rand_index), quote = FALSE)
  print(paste("Number of iterations:", biclust_result$n_iterations), quote = FALSE)
  print(paste("Silhoutte of first disturbed cluster likelihood (aka how complex was the first likelihood):", biclust_result$db), quote = FALSE)
  print(paste("BIC_all:", biclust_result$BIC), quote = FALSE)
  print(paste("taget_genes_residual_var:"), quote = FALSE)
  print(biclust_result$taget_genes_residual_var, quote = FALSE)
}
