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
# source(file.path(R_path, "plot_cluster_history.R"))
# source(file.path(R_path, "plot_loglikelihood.R"))


# simpler regression function:
weighted_lm <- function(X, Y, W = diag(nrow(X))) {

  check <- function(X) {
    if (!is.matrix(X)) {
      # Try to coerce X into a matrix
      X <- try(as.matrix(X), silent = TRUE)
      if (!is.matrix(X)) {
        stop("Input not coercible to matrix")
      }
    }
    return(X)
  }

  # make sure X,Y,W either are or can be coerced to matrices and if not return error
  check(X); check(Y); check(W)
  BETAS <- pracma::pinv(t(X) %*% W %*% X) %*%
    t(X) %*%
    W %*%
    Y

  #returns betas, one column per target variable
  return(BETAS)
}

ridge_weighted_lm_ <- function(X, Y, W = diag(nrow(X)), lambda = 1 / nrow(X)) {
  # Estimates betas using weighted rigde regression, requires additional parameter lambda

  BETAS <- weighted_lm(X, Y, W)
  N <- nrow(X)
  BETAS_ridge <- BETAS / (1 + N * lambda)

  return(BETAS_ridge)
}

test_weighted_lm <- function(X, Y) {

  lm_result <- lm(Y ~ 0 + X)$coefficients #intercept-free method
  my_result <- weighted_lm(X, Y, W = diag(nrow(X)))

  mean(abs(unlist(lm_result) - my_result))
}

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
biclust_lm <- function(dat = dat,
                    dat_test = NULL,
                    max_iter = 50,
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


  # TODO: Split into train and test here instead of in script
  n_total_cells_test    <- nrow(dat_test)

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
  weights_all <- vector(mode = "list", length = max_iter)
  BIC_all <- vector(mode = "list", length = max_iter)
  target_genes_residual_var_all <- vector(mode = "list", length = max_iter)
  # Set the current cell clustering
  current_cell_cluster_allocation <- initial_clustering

  stop_iterating_flag <- 0  # Flag if we have converged
  for (i_main in 1:max_iter) {
    start.time <- Sys.time()

    ################################################################
    ##### M-step, compute estimates for \pi_k and model parameters #
    ################################################################

    ###### Fit model to each cell cluster ####
    ###### M.1                           ####
    models <- vector(mode = "list", length = n_cell_clusters)

    # #dev
    # models2 <- vector(mode = "list", length = n_cell_clusters)
    # models_eval <- vector(mode = "list", length = n_cell_clusters)

    for (i_cell_cluster in 1:n_cell_clusters) {
      cell_cluster_rows <- which(current_cell_cluster_allocation == i_cell_cluster)
      cell_cluster_target_genes <- as.matrix(dat[cell_cluster_rows, ind_targetgenes, drop = FALSE])
      cell_cluster_regulator_genes <- as.matrix(dat[cell_cluster_rows, ind_reggenes, drop = FALSE])

      if (nrow(cell_cluster_target_genes) == 0) {
        print(paste("Number of cells in cell cluster ", i_cell_cluster, "is 0"), quote = FALSE)
        return(NA)
      }

      if (use_weights == FALSE || i_main == 1) {
        current_weights <- diag(nrow(cell_cluster_regulator_genes))
      }
      else if (use_weights == TRUE && i_main == 1) {
        current_weights <- diag(nrow(cell_cluster_regulator_genes)) * cluster_proportions[i_cell_cluster]
      }
      else {
        # Current weights is a n_cell x n_cell matrix with the weigths on the diagonal
        current_weights <- weights_all[[i_main - 1]][cell_cluster_rows, drop = FALSE]

        # current_weights <- current_weights / sum(current_weights)
        if (length(current_weights) == 1) {
          current_weights <- matrix(current_weights, nrow = 1, ncol = 1)
        }else if (length(current_weights) == 0) {
          print("Current weights are:", quote = FALSE)
          print(current_weights, quote = FALSE)
          stop("Length of weights is 0")

        }else {
          current_weights <- diag(current_weights)
        }

      }
      # models2[[i_cell_cluster]] <- lm(formula = 'cell_cluster_target_genes ~ 0 + cell_cluster_regulator_genes',
      #                                 data = environment())$coefficients

      # Here we use a ridge weighted LM, BUT IT COULD BE ANYTHING,
      #    the only important part is that it allows for assignment of cluster probabilities to observations.
      # This returns the parameters as a matrix
      models[[i_cell_cluster]] <- ridge_weighted_lm_(
        X = cell_cluster_regulator_genes,
        Y = cell_cluster_target_genes,
        W = current_weights,
        lambda = penalization_lambda
      )

      #  cat(
      #    mean(dat[cell_cluster_rows,2] != as.numeric(initial_clustering[cell_cluster_rows])), '\n'
      #  )


      # models_eval[[i_cell_cluster]]  <- test_weighted_lm(
      #   X = as.matrix(dat[cell_cluster_rows, ind_reggenes]),
      #   Y = as.matrix(dat[cell_cluster_rows, ind_targetgenes]))

    }


    # dev
    # target_genes_residual_var2 <- matrix(data = 0, nrow = n_cell_clusters, ncol = n_target_genes)
    # dat is dat <- cbind(target_expression, regulator_expression), e.g. a 2x100, with e.g. the first 50 rows being true cell cluster 1
    # 100x2 * 2x1

    #### calculate one variance for each target gene given each model ####
    ###### M.2                                                        ####

    # Pre-allocate residual variance estimates
    target_genes_residual_var <- matrix(data = 0, nrow = n_cell_clusters, ncol = n_target_genes)

    # Output: n_cell_clusters x n_target_genes
    for (i_cell_cluster in seq_len(n_cell_clusters)) {
      current_rows <- which(current_cell_cluster_allocation == i_cell_cluster)
      current_regulator_genes <- as.matrix(dat[current_rows, ind_reggenes, drop = FALSE])
      current_target_genes <- as.matrix(dat[current_rows, ind_targetgenes, drop = FALSE])
      cell_cluster_betas <- models[[i_cell_cluster]] # $coefficients # with our own lm we only save the coeffs anyway


      predicted_values <- current_regulator_genes %*% cell_cluster_betas

      residuals <- current_target_genes - predicted_values
      # target_genes_residual_var[i_cell_cluster,] <- diag(var(residuals))  # maybe not necessary to calculate entire matrix
      target_genes_residual_var[i_cell_cluster,] <- colSums(residuals^2) / (length(current_rows) - 1)

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

    weights <- sweep(exp(likelihood), 2, cluster_proportions, "*")


    # weights_2 <- sweep((likelihood), 2, cluster_proportions, "*")
    # weights_2 <- sweep(weights_2, 1, rowSums(weights_2), "/")
    likelihood <- weights

    #
    # weights_2 <- sweep(likelihood, 2, cluster_proportions, "+")
    # weights_2 <- sweep(weights_2, 2, colSums(weights_2), "-")
    # # print(str(weights_2))
    # # print(str(weights))
    #
    weights_all[[i_main]] <- weights

    likelihood_all[[i_main]] <- likelihood



    ###### compute BIC
    #first, predict cluster allocation
    likelihood_test <- loglikelihood_calc_matrix(dat_test,
                                                 models,
                                                 target_genes_residual_var,
                                                 n_cell_clusters,
                                                 ind_reggenes,
                                                 ind_targetgenes)

    predicted_cluster_allocation_test <- sapply(seq_len(nrow(likelihood_test)),
                                                function(row) which.max(likelihood_test[row,]))

    # compute cluster allocation proportions
    cluster_proportions_test <- vector(length = n_cell_clusters)
    for (i_cell_cluster in seq_len(n_cell_clusters)) {
      current_cluster_proportion <- sum(predicted_cluster_allocation_test == i_cell_cluster) / length(predicted_cluster_allocation_test)
      if (current_cluster_proportion == 0) {
        cluster_proportions_test[i_cell_cluster] <- 0.0001
        warning("Test data has zero cells in predicted cluster")
      }else
        cluster_proportions_test[i_cell_cluster] <- current_cluster_proportion
    }
    #calculate likelihood of model for test data

    weights_test <- sweep(exp(likelihood_test), 2, cluster_proportions_test, "*")


    big_logLhat_in_BIC <- sum(log(rowSums(weights_test)))
    k_in_BIC <- (n_target_genes + n_regulator_genes) * n_cell_clusters
    n_in_BIC <- n_total_cells_test

    BIC <- k_in_BIC * log(n_in_BIC) - 2 * big_logLhat_in_BIC
    BIC_all[[i_main]] <- BIC

    if (i_main == 1) {
      # no_factor_cluster <- as.numeric(levels(initial_clustering))[initial_clustering]
      # db <- index.DB(likelihood, no_factor_cluster)$DB
      # db <- mean(silhouette(no_factor_cluster, dist(likelihood))[, 3])
      db <- 0 # hack
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

    if(length(unique(updated_cell_clust)) == 1){
      print("Everything is put into the same cell cluster")
    }
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

    rand_index_true_cluster <- aricode::RI(dat$true_cell_cluster_allocation, updated_cell_clust)
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
                               true_cell_cluster_allocation_vector = dat$true_cell_cluster_allocation)
    hist_plot_loglikelihood(dat,
                            likelihood,
                            n_cell_clusters,
                            penalization_lambda,
                            output_path,
                            i_main)
  }

  start.time <- Sys.time()

  cell_cluster_history_plotting <- cbind('Cell ID' = cell_cluster_history[, 1],
                                         'True cell cluster allocation' = dat$true_cell_cluster_allocation,
                                         cell_cluster_history[, 2:ncol(cell_cluster_history)])
  png(file.path(output_path, paste0("lm_Alluvial_diagram_lambda_", round(penalization_lambda, 3), ".png")),
      width = 1024 + ncol(cell_cluster_history_plotting) * 40, height = 1024, units = "px", res = 150)
  plot_cluster_history(cell_cluster_history = cell_cluster_history_plotting, correct_plot = FALSE)
  dev.off()
  time_taken <- round(Sys.time() - start.time, 2)
  print(paste(" Iterations complete, Alluvial plot took", time_taken, "seconds."), quote = FALSE)

  return(list("rand_index" = rand_index_true_cluster,
              "n_iterations" = i_main,
              "db" = db,
              "BIC" = BIC_all,
              "target_genes_residual_var" = target_genes_residual_var_all))
}


# Runs only when script is run by itself
# || interactive()
if (sys.nframe() == 0) {
  # ... do main stuff

  # Set seed for example
  set.seed(1234)


  #############################################
  ############ data for dev ###################
  #############################################

  R_path <- here::here("R")
  # source(file.path(R_path, "generate_data_lm.R"))
  # source(file.path(R_path, "randomise_cluster_labels.R"))


  dat <- generate_data_lm(n_cell_clusters = 3,
                          n_target_gene_type = 2,  # We have x named target genes that have one expression per cell
                          n_regulator_gene_type = 2,  # We have x named regulator genes that have one expression per cell
                          n_cells = c(1000, 5000, 10000),
                          regulator_means = c(1, 2, 5),  # Regulator mean expression in each cell cluster.
                          regulator_standard_deviations = c(0.1, 0.2, 0.3),  # Regulator sd for expression in each cell cluster.
                          coefficients_standard_deviation = 100, # 'The betas/slopes'. One per target gene. Instead of providing mean and standard deviation for each target gene, we provide the standard deviation from which these will be generated. Mean will be 0.
                          target_gene_type_standard_deviation = 3
  )

  true_betas <- dat$true_betas
  dat <- dat$dat
  # dat[, 'true_cell_cluster_allocation'] <- paste("Cluster", pull(dat, 'true_cell_cluster_allocation'))
  # These needs to be strings for discrete labels in pca plot
  # that fucks up the code, though
  pca_res <- prcomp(dat[, 3:ncol(dat)], scale. = TRUE)
  p <- ggplot2::autoplot(pca_res, data = dat, colour = 'true_cell_cluster_allocation')


  n_total_cells <- nrow(dat)
  ind_targetgenes <- which(str_detect(colnames(dat), "t\\d"))
  ind_reggenes <- which(str_detect(colnames(dat), "r\\d"))

  disturbed_initial_cell_clust <- randomise_cluster_labels(cluster_labels = dat$true_cell_cluster_allocation,
                                                           fraction_randomised = 10 / 100)

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
  penalization_lambda <- 0.000000001
  i_cell_cluster <- 1
  i_main <- 1

  use_weights = TRUE
  use_complex_cluster_allocation = FALSE

  demo_path <- here::here("demo")
  output_path <- demo_path


  ###########################################
  ###END initialise variables for dev #######
  ###########################################
  biclust(dat = dat,
          max_iter = 50,
          initial_clustering,
          n_target_genes,
          n_regulator_genes,
          n_total_cells,
          n_cell_clusters,
          ind_targetgenes,
          ind_reggenes,
          output_path,
          penalization_lambda = 0.5)

}
