#!/usr/bin/Rscript

library(tidyverse)
library(aricode)     # To calculate rand index
library(reshape2)
library(here)        # To work with paths
library(pracma)      # For pseudo inverse
library(stats)       # For kmeans
library(ppclust)
library(clusterSim)  # for db
library(Rmpfr)  # for small loglikihoods
library(parallel)
library(doParallel)
library(foreach)

# Get absolute path where script is located, by using relative paths.
R_path <- here::here("R")
source(file.path(R_path, "plot_cluster_history.R"))
source(file.path(R_path, "plot_loglikelihood.R"))

hms_span <- function(start, end) {
  dsec <- as.numeric(difftime(end, start, unit = "secs"))
  hours <- floor(dsec / 3600)
  minutes <- floor((dsec - 3600 * hours) / 60)
  seconds <- dsec - 3600 * hours - 60 * minutes
  paste0(
    sapply(c(hours, minutes, seconds), function(x) {
      formatC(x, width = 2, format = "d", flag = "0")
    }), collapse = ":")
}

# Because scregclust doesn't want a index vector but a logical one
# and also in loglike calc we need the inverse of a index vector
inverse_which <- function(indices, output_length)
{
  rv <- logical(output_length)
  if (length(indices) > 0)
  {
    rv[indices] <- TRUE
  }
  return(rv)
}


loglikelihood_calc_matrix <- function(all_genes,
                                      models,
                                      target_genes_residual_var,
                                      n_cell_clusters,
                                      ind_reggenes,
                                      ind_targetgenes,
                                      data_split_for_scregclust,
                                      current_cell_cluster_allocation) {
  # For all cells, calculate the likelihood of coming from the model corresponding to each
  n_cells <- nrow(all_genes)
  loglikelihood <- matrix(data = 0, nrow = n_cells, ncol = n_cell_clusters)

  regulator_genes <- as.matrix(all_genes[, ind_reggenes, drop = FALSE])  # 1xr -> cxr
  target_genes <- as.matrix(all_genes[, ind_targetgenes, drop = FALSE])  # 1xt -> cxt

  parameters <- vector(mode = "list", length = n_cell_clusters)

  for (i_cell_cluster in seq_len(n_cell_clusters)) {
    # Standardize data like in scregclust
    cell_cluster_rows_ind <- which(current_cell_cluster_allocation == i_cell_cluster)

    training_data_ind <- which(data_split_for_scregclust[[i_cell_cluster]] == 1)

    training_data_for_scregclust_ind <- cell_cluster_rows_ind[training_data_ind]

    rest_of_data_ind <- which(!inverse_which(training_data_for_scregclust_ind, n_cells))

    #save this function and return it for use in the predict function

    res <- standardize_like_scregclust(xvals = regulator_genes, yvals = target_genes,
                                       training_data_ind = training_data_for_scregclust_ind,
                                       test_data_ind = rest_of_data_ind)

    parameters[[i_cell_cluster]] <- res$params

    regulator_genes_normalized <- res[['xvals']]

    target_genes_normalized <- res[['yvals']]

    #--------------------

    print(paste("  Calculating loglikelihood for cell cluster", i_cell_cluster), quote = FALSE)


    cell_cluster_betas <- models[[i_cell_cluster]] # rxt

    predicted_value <- regulator_genes_normalized %*% cell_cluster_betas  # 1xt -> cxt
    cell_squared_error <- (target_genes_normalized - predicted_value)^2  # 1xt -> cxt
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
  return(list(
    'loglikelihood' = loglikelihood,
    'parameters' = parameters
  )
  )
}


standardize_like_scregclust <- function(xvals, yvals, training_data_ind, test_data_ind) {


  org_yvals <- yvals
  org_xvals <- xvals

  yvals[training_data_ind,] <- scale(yvals[training_data_ind,], scale = FALSE)
  xvals[training_data_ind,] <- scale(xvals[training_data_ind,])

  yval_colmeans  <- colMeans(org_yvals[training_data_ind,])
  xvals_colmeans <- colMeans(org_xvals[training_data_ind,])
  xval_scale     <- apply(org_xvals[training_data_ind,], 2, sd)

  yvals[test_data_ind,] <- scale(yvals[test_data_ind,], yval_colmeans, scale = FALSE)
  xvals[test_data_ind,] <- scale(xvals[test_data_ind,], xvals_colmeans, xval_scale )

  return(list("xvals" = xvals, "yvals" = yvals,
              "params" = list(
                'yval_colmeans'  = yval_colmeans,
                'xvals_colmeans' = xvals_colmeans,
                'xval_scale'     = xval_scale
              )
  )
  )
}


#' Bisc
#'
#' Cluster cells based with Classification EM based on lm with penality
#' @param dat Main data
#' @param cell_id
#' @param true_cell_cluster_allocation, The tru cell cluster allocation for Rand index caluclations
#' @param max_iter Max number of iterations to run
#' @param n_target_gene_clusters
#' @param initial_clustering The start clustering
#' @param n_cell_clusters Number of cell clusters
#' @param ind_targetgenes Indexes of where the target gene columns in dat are
#' @param ind_reggenes
#' @param output_path Output path of plots
#' @param penalization_lambda The penalization lambda
#' @param calculate_optimization_target [TRUE] If to calculate and return eq 2.2 and 2.4 in the CEM paper. 2.4 is basically what the algo is optimizing.
#' @param calculate_silhoutte [FALSE]
#' @param calculate_davies_bouldin_index [FALSE]
#' @param plot_suffix [""]
#' @param always_use_flat_prior [FALSE]
#' @param use_garbage_cluster_targets [FALSE]
#' @param retain_gene_clusters [TRUE] If to remember the gene module allocation from last iteration as to try and speed up scregclust.
#' @return Nothing yet, maybe cluster labels
#' @export
#'
bisc <- function(
    dat = dat,
    cell_id,
    true_cell_cluster_allocation,
    max_iter = 50,
    n_target_gene_clusters,
    initial_clustering,
    n_cell_clusters,
    ind_targetgenes,
    ind_reggenes,
    output_path,
    penalization_lambda = 0.05,
    calculate_optimization_target = TRUE,
    calculate_silhoutte = FALSE,
    calculate_davies_bouldin_index = FALSE,
    plot_suffix = "",
    always_use_flat_prior = FALSE,
    use_garbage_cluster_targets = FALSE,
    retain_gene_clusters = TRUE
) {

  if (!is.factor(initial_clustering)) {
    stop("The variable initial_clustering needs to be a factor vector.")
  }
  n_target_genes <- length(ind_targetgenes)
  n_regulator_genes <- length(ind_reggenes)
  converged <- FALSE
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
  RIs_all <- vector(mode = "numeric", length = max_iter)
  max_loglikelihood_all <- vector(mode = "numeric", length = max_iter)
  CML_C2_all <- vector(mode = "numeric", length = max_iter)
  target_genes_residual_var_all <- vector(mode = "list", length = max_iter)
  # Set the current cell clustering
  current_cell_cluster_allocation <- initial_clustering

  # preallocate gene cluster tracker
  previous_gene_modules = vector(mode = "list", length = n_cell_clusters)

  stop_iterating_flag <- 0  # Flag if we have converged
  # (results1[[1]][[3]]$results[[1]]$output[[1]]$models)
  scregclust_final_result_models <- vector(mode = "list", length = n_cell_clusters)
  # (results1[[1]][[3]]$results[[1]]$output[[1]]$module)
  scregclust_final_result_module <- vector(mode = "list", length = n_cell_clusters)
  scregclust_final_result_coeffs <- vector(mode = "list", length = n_cell_clusters)
  scregclust_final_result_signs <- vector(mode = "list", length = n_cell_clusters)
  scregclust_final_result_reg_table <- vector(mode = "list", length = n_cell_clusters)


  for (i_main in 1:max_iter) {
    print(paste(" Running iteration", i_main), quote = FALSE)
    start_time_iteration <- Sys.time()

    ################################################################
    ##### M-step, compute estimates for \pi_k and model parameters #
    ##### \pi_k are cluster proportions                            #
    ##### model parameters are betas and signas from SCREG         #
    ################################################################

    ###### Fit model to each cell cluster ####
    ###### M.1                           ####
    models <- vector(mode = "list", length = n_cell_clusters)
    data_split_for_scregclust <- vector(mode = "list", length = n_cell_clusters)

    for (i_cell_cluster in 1:n_cell_clusters) {
      cell_cluster_rows <- which(current_cell_cluster_allocation == i_cell_cluster)
      cell_cluster_target_genes <- as.matrix(dat[cell_cluster_rows, ind_targetgenes, drop = FALSE])
      cell_cluster_regulator_genes <- as.matrix(dat[cell_cluster_rows, ind_reggenes, drop = FALSE])


      print(paste("  Calculating betas for cell cluster", i_cell_cluster, "with", nrow(cell_cluster_target_genes), "cells."), quote = FALSE)

      if (nrow(cell_cluster_target_genes) == 0) {
        print(paste("Number of cells in cell cluster ", i_cell_cluster, "is 0"), quote = FALSE)
        return(NA)
      }

      #  scregclust is function that takes expression = p rows of genes and n columns of cells.
      screg_out_betas <- matrix(data = 0, nrow = n_regulator_genes, ncol = n_target_genes)

      if (!nrow(cell_cluster_target_genes) == 0) {
        indata_for_scregclust <- rbind(t(cell_cluster_target_genes), t(cell_cluster_regulator_genes))

        # Training data are represented by 1 and test data by 2
        data_split_for_scregclust[[i_cell_cluster]] <- sample(rep(1:2, length.out = ncol(indata_for_scregclust)))

        n_training_data <- sum(data_split_for_scregclust[[i_cell_cluster]] == 1)
        n_test_data <- sum(data_split_for_scregclust[[i_cell_cluster]] == 2)

        # if ((length(ind_reggenes) > n_test_data) || (length(ind_reggenes) > n_training_data)) {
        #   print(paste("   The split into training/test data for scregclust resulted in less cells than there are regulators."), quote = FALSE)
        #   print(paste("    There are", length(ind_reggenes), "regulator genes."), quote = FALSE)
        #   print(paste0("    There are ", ncol(indata_for_scregclust), " cells in cell cluster ", i_cell_cluster, "."), quote = FALSE)
        #   print(paste("    There are", n_test_data, "cells in the test data."), quote = FALSE)
        #   print(paste("    There are", n_training_data, "cells in the training data."), quote = FALSE)
        #   return(NA)
        # }

        # Fix constant genes
        ind_cells_train <- which(data_split_for_scregclust[[i_cell_cluster]] == 1)
        ind_cells_test <- which(data_split_for_scregclust[[i_cell_cluster]] == 2)
        # print(str(indata_for_scregclust))

        constant_ind_test <- which(apply(indata_for_scregclust[, ind_cells_test, drop = FALSE], 1, sd) == 0)

        non_constant_genes_ind_test <- which(apply(indata_for_scregclust[, ind_cells_test, drop = FALSE], 1, sd) != 0)
        non_constant_genes_ind_train <- which(apply(indata_for_scregclust[, ind_cells_train, drop = FALSE], 1, sd) != 0)

        non_constant_ind_genes <- sort(intersect(non_constant_genes_ind_test, non_constant_genes_ind_train))
        non_constant_ind_regulator_genes <- non_constant_ind_genes[non_constant_ind_genes > n_target_genes] - n_target_genes


        if(retain_gene_clusters){
          if(i_main == 1){

            # sink(tempfile()) # Shut scregclust up
            screg_out <- scregclust::scregclust(
              expression = indata_for_scregclust,  # p rows of genes, n columns of cells
              split_indices = data_split_for_scregclust[[i_cell_cluster]],
              genesymbols = 1:(n_target_genes + n_regulator_genes),  # Gene row numbers
              is_regulator = inverse_which(indices = ind_reggenes, output_length = n_regulator_genes + n_target_genes),  # Vector indicating which genes are regulators
              n_modules = n_target_gene_clusters[[i_cell_cluster]],
              penalization = penalization_lambda,
              verbose = FALSE,
              max_optim_iter = 100000L,
              # tol_coop_rel = 1e-06,
              # tol_coop_abs = 1e-06,
              # tol_nnls = 1e-02,
              min_module_size = 1L,
              n_cycles = 25,
              center = FALSE,
            )
            # sink()

            if(use_garbage_cluster_targets){
              previous_gene_modules[[i_cell_cluster]] <- screg_out$results[[1]]$output[[1]]$module
            }else{
              previous_gene_modules[[i_cell_cluster]] <- screg_out$results[[1]]$output[[1]]$module_all
            }

          }else{

            is_regulator <- inverse_which(indices = ind_reggenes, output_length = n_regulator_genes + n_target_genes)

            # sink(tempfile()) # Shut scregclust up
            screg_out <- scregclust::scregclust(
              expression = indata_for_scregclust,  # p rows of genes, n columns of cells
              split_indices = data_split_for_scregclust[[i_cell_cluster]],
              genesymbols = 1:(n_target_genes + n_regulator_genes),  # Gene row numbers
              is_regulator = is_regulator,  # Vector indicating which genes are regulators
              n_modules = n_target_gene_clusters[[i_cell_cluster]],
              penalization = penalization_lambda,
              initial_target_modules = previous_gene_modules[[i_cell_cluster]][!is_regulator],
              verbose = FALSE,
              n_cycles = 25,
              center = FALSE,
              max_optim_iter = 100000L,
              # tol_coop_rel = 1e-06,
              # tol_coop_abs = 1e-06,
              # tol_nnls = 1e-02,
              min_module_size = 1L
            )
            # sink()

          }

        }else{
          screg_out <- scregclust::scregclust(
            expression = indata_for_scregclust,  # p rows of genes, n columns of cells
            split_indices = data_split_for_scregclust[[i_cell_cluster]],
            genesymbols = 1:(n_target_genes + n_regulator_genes),  # Gene row numbers
            is_regulator = inverse_which(indices = ind_reggenes, output_length = n_regulator_genes + n_target_genes),  # Vector indicating which genes are regulators
            n_modules = n_target_gene_clusters[[i_cell_cluster]],
            penalization = penalization_lambda,
            verbose = FALSE,
            n_cycles = 25,
            center = FALSE,
            max_optim_iter = 100000L,
            # tol_coop_rel = 1e-06,
            # tol_coop_abs = 1e-06,
            # tol_nnls = 1e-02,
            min_module_size = 1L
          )
        }

        if(use_garbage_cluster_targets){
          scregclust_final_result_module[[i_cell_cluster]] <- screg_out$results[[1]]$output[[1]]$module
        }else{
          scregclust_final_result_module[[i_cell_cluster]] <- screg_out$results[[1]]$output[[1]]$module_all
        }
        scregclust_final_result_models[[i_cell_cluster]] <- screg_out$results[[1]]$output[[1]]$models
        scregclust_final_result_coeffs[[i_cell_cluster]] <- screg_out$results[[1]]$output[[1]]$coeffs
        scregclust_final_result_signs[[i_cell_cluster]] <- screg_out$results[[1]]$output[[1]]$signs
        scregclust_final_result_reg_table[[i_cell_cluster]] <- screg_out$results[[1]]$output[[1]]$reg_table


        # v Step 1: Regulators selected (in 8ms)
        # # Final counts
        #      Cluster  Noise |     1 |     2
        #      Targets   1627 |     0 |     0
        #   Regulators      - |     0 |     0
        #
        # v Step 2a: Coefficients re-estimated (in 8ms)
        # v Post-processing and packaging done (in 4ms)
        #
        # v Finished clustering with penalization = 0.5
        #
        # v Finished!
        # Warning in rm(beta_logical, coeffs_new_style, temp_coeffs, current_beta_logical) :
        #   object 'temp_coeffs' not found
        # Warning in rm(beta_logical, coeffs_new_style, temp_coeffs, current_beta_logical) :
        #   object 'current_beta_logical' not found
        # [1] Cell cluster 1 betas is NULL for scregclust output.
        # [1]
        # [1]
        # [1] penalization_lambda 0.5 is NA # TODO

        # Since the output in scregclust changed with after commit e0516fcdd89a6a549d7906f6502f41586c3ed47f
        # we need to convert the current output to the old output style
        print(paste("  Converting scregclust output to old style."), quote = FALSE)
        beta_logical <- screg_out$results[[1]]$output[[1]]$models
        coeffs_old_style <- vector(mode = "list", length = ncol(beta_logical))
        coeffs_new_style <- screg_out$results[[1]]$output[[1]]$coeffs
        for (i_target_gene_cluster in seq(length(coeffs_old_style))) {
          if (is.null(coeffs_new_style[[i_target_gene_cluster]])) {
            coeffs_old_style[[i_target_gene_cluster]] <- NULL
            temp_coeffs <- NULL
            current_beta_logical <- NULL
          }else {
            temp_coeffs <- matrix(0, nrow = nrow(beta_logical), ncol = ncol(coeffs_new_style[[i_target_gene_cluster]]))
            current_beta_logical <- beta_logical[, i_target_gene_cluster]  # Reduces to a vector in R
            temp_coeffs[current_beta_logical,] <- coeffs_new_style[[i_target_gene_cluster]]
            coeffs_old_style[[i_target_gene_cluster]] <- temp_coeffs
          }
        }
        rm(beta_logical, coeffs_new_style, temp_coeffs, current_beta_logical)


        # This is OK, things work out anyway
        # if (!screg_out$results[[1]]$converged) {
        #   print(paste("scregclust didn't converge for cell cluster", i_cell_cluster,"."), quote = FALSE)
        #   return(NA)
        # }

        screg_out_betas_temp <- do.call(cbind, coeffs_old_style)  # Merge betas into one matrix
        if (is.null(screg_out_betas_temp)) {
          print(paste("Cell cluster", i_cell_cluster, "betas is NULL for scregclust output."), quote = FALSE)
          return(NA)
        }

        screg_out_betas <- matrix(NA, nrow = n_regulator_genes, ncol = ncol(screg_out_betas_temp))
        screg_out_betas[non_constant_ind_regulator_genes,] <- screg_out_betas_temp


        # if (length(screg_out$results[[1]]$output[[1]]$cluster) < (n_target_genes + n_regulator_genes)) {
        #   print(paste("Number of outputted genes with named clusters is less than we fed into scregclust.",
        #               length(screg_out$results[[1]]$output[[1]]$cluster),
        #               "is less than target + reg =",
        #               n_target_genes, "+", n_regulator_genes, "=", n_regulator_genes + n_target_genes), quote = FALSE)
        #   return(NA)
        # }

        # Create cluster_indexes which maps columns in screg_out_betas to correct places under some assumptions of order:
        # target_gene_cluster_names, e.g. 2 2 1 3 1 3,
        # becomes cluster_indexes: 3 4 1 5 2 6
        # CAN ALSO BE DONE WITH THE FOLLOWING CODE THAT WE HAVE USED BEFORE:
        # indexes_of_not_deleted_target_gene_clusters <- which(sapply(out_list[[i_target_cell_cluster]]$results[[1]]$output[[1]]$coeffs, FUN=function(x) !is.null(x)))
        # target_gene_cluster_index <- out_list[[i_target_cell_cluster]]$results[[1]]$output[[1]]$cluster[1:n_target_genes]  # n_target_genes, e.g. 1 1 1 1 3 3 3 3 3 3 3 3 1 3 1 3 3 3 3 3 1 3 1 1 3 1 3 3 1 3
        # target_gene_cluster_index <- unlist(sapply(indexes_of_not_deleted_target_gene_clusters, FUN=function(x) which(target_gene_cluster_index==x)))

        # This is a work around because sregclust removes constant genes before it starts
        target_gene_cluster_names <- rep(-1, n_target_genes + n_regulator_genes)
        if (use_garbage_cluster_targets) {
          target_gene_cluster_names_temp <- screg_out$results[[1]]$output[[1]]$module  # [1:n_target_genes]
        }else {
          target_gene_cluster_names_temp <- screg_out$results[[1]]$output[[1]]$module_all  # [1:n_target_genes]
        }
        target_gene_cluster_names[non_constant_ind_genes] <- target_gene_cluster_names_temp
        target_gene_cluster_names <- target_gene_cluster_names[1:n_target_genes]

        # Fixa den här. target_gene_cluster_names är 47. Ut kommer cluster_indexes som är 47
        cluster_indexes <- rep(NA, length(target_gene_cluster_names))
        tmp_max_ind <- 0
        for (i_target_gene_cluster in 1:n_target_gene_clusters[[i_cell_cluster]]) {
          ind_tmp <- which(target_gene_cluster_names == i_target_gene_cluster)
          n_target_genes_in_cluster <- length(ind_tmp)
          cluster_indexes[ind_tmp] <- 1:n_target_genes_in_cluster + tmp_max_ind  # Get names of clusters for each target gene
          tmp_max_ind <- n_target_genes_in_cluster
        }
        # print(max(cluster_indexes, na.rm = TRUE), quote = FALSE)
        # print(ncol(screg_out_betas), quote = FALSE)

        screg_out_betas <- screg_out_betas[, cluster_indexes]  # Sort the matrix


      }
      # screg_out_sigmas <- do.call(c, screg_out$results[[1]]$output[[1]]$sigmas)
      # screg_out_sigmas <- screg_out_sigmas[cluster_indexes]

      models[[i_cell_cluster]] <- screg_out_betas


      # TODO: alternatively write some predict function and use extracted variance estimates
    }

    # scregclust_final_result_coeffs_neat <- models

    #### calculate one variance for each target gene given each model ####
    ###### M.2                                                        ####

    # Pre-allocate residual variance estimates
    target_genes_residual_var <- matrix(data = 0, nrow = n_cell_clusters, ncol = n_target_genes)

    # Output: n_cell_clusters x n_target_genes
    for (i_cell_cluster in seq_len(n_cell_clusters)) {
      current_rows <- which(current_cell_cluster_allocation == i_cell_cluster)
      current_regulator_genes <- as.matrix(dat[current_rows, ind_reggenes, drop = FALSE])
      current_target_genes <- as.matrix(dat[current_rows, ind_targetgenes, drop = FALSE])
      print(paste("  Calculating residual variance for cell cluster", i_cell_cluster, "with", length(current_rows), "cells."), quote = FALSE)

      # Standardize
      training_data_ind <- which(data_split_for_scregclust[[i_cell_cluster]] == 1)
      test_data_ind <- which(data_split_for_scregclust[[i_cell_cluster]] == 2)
      res <- standardize_like_scregclust(xvals = current_regulator_genes, yvals = current_target_genes,
                                         training_data_ind = training_data_ind, test_data_ind = test_data_ind)
      current_regulator_genes <- res[['xvals']]
      current_target_genes <- res[['yvals']]


      cell_cluster_betas <- models[[i_cell_cluster]] # $coefficients # with our own lm we only save the coeffs anyway

      predicted_values <- current_regulator_genes %*% cell_cluster_betas

      residuals <- current_target_genes - predicted_values

      # target_genes_residual_var[i_cell_cluster,] <- diag(var(residuals))  # maybe not necessary to calculate entire matrix
      # TODO: this should be possible to just extract from SCREG, do that and compare to make sure it's approx same
      tmp_res_var <- colSums(residuals^2) / (length(current_rows) - 1)
      tmp_res_var[is.na(tmp_res_var)] <- mean(tmp_res_var, na.rm = TRUE)  # Make target genes that were thrown in the rubbish cluster have large variance
      target_genes_residual_var[i_cell_cluster,] <- tmp_res_var

    }
    target_genes_residual_var_all[[i_main]] <- target_genes_residual_var


    # Calculated loglikelihoods
    # M.3

    loglikelihood <- loglikelihood_calc_matrix(all_genes = dat,
                                               models = models,
                                               target_genes_residual_var = target_genes_residual_var,
                                               n_cell_clusters = n_cell_clusters,
                                               ind_reggenes = ind_reggenes,
                                               ind_targetgenes = ind_targetgenes,
                                               data_split_for_scregclust = data_split_for_scregclust,
                                               current_cell_cluster_allocation = current_cell_cluster_allocation)


    # separate out normalisation parameters
    normalisation_parameters <- loglikelihood$parameters
    loglikelihood            <- loglikelihood$loglikelihood


    if (!always_use_flat_prior) {
      cluster_proportions <- vector(length = n_cell_clusters)
      for (i_cell_cluster in seq_len(n_cell_clusters)) {
        print(paste("  Calculating cluster proportions for cell cluster", i_cell_cluster), quote = FALSE)
        current_cluster_proportion <- sum(current_cell_cluster_allocation == i_cell_cluster) / length(current_cell_cluster_allocation)

        # If some cluster has been completely omitted, give it a nonzero proportion anyway for next iteration
        # Even without posterior_probability this is good to include for cluster allocation (though it is arbitrary)
        if (current_cluster_proportion == 0) {
          cluster_proportions[i_cell_cluster] <- 0.0001
          stop_iterating_flag <- T
        }else
          cluster_proportions[i_cell_cluster] <- current_cluster_proportion
      }
    } else {
      cluster_proportions <- (vector(length = n_cell_clusters) + 1) / n_cell_clusters
    }


    ####################################################################
    ##### E-step #######################################################
    ####### Compute posterior probabilities ############################
    ####### For each cell to belong to each cluster ####################

    # If everything goes well this should work the same for screg as with our
    # lm examples


    print(paste("  Doing final posterior probability calculations"), quote = FALSE)
    loglikelihood <- Rmpfr::mpfr(x = loglikelihood, precBits = 512)  # cells x cell_clusters

    # unnormalized_joint_posterior is a matrix (cells x cell_clusters). cluster_proportions is a vector cell_clusters long
    unnormalized_joint_posterior <- sweep(exp(loglikelihood), 2, (cluster_proportions), "*")

    # --------------------
    # Normalize posterior likelihood
    # Below parallise the this calculation (e.g speeds up from 1:32 min to 12 sec)
    # posterior_probability <- sweep(posterior_probability, 1, rowSums(posterior_probability), "/")

    # Define the number of cores to use
    num_cores <- max(detectCores() - 2, 1)

    # # Create a parallel backend
    cl <- parallel::makeCluster(num_cores,  type = "PSOCK")

    # Register the parallel backend
    registerDoParallel(cl)
    num_rows <- nrow(unnormalized_joint_posterior)

    split_indices <- split(1:num_rows, cut(1:num_rows, num_cores))

    # Define the function to parallelize
    parallel_sweep <- function(x) {
      Rmpfr::asNumeric(sweep(x, 1, rowSums(x), "/"))
    }

    # Parallelize the operation using foreach
    normalized_posterior_probability <- foreach(i = split_indices, .packages = c("Rmpfr"), .combine = 'rbind') %dopar% {
      parallel_sweep(unnormalized_joint_posterior[i, , drop = FALSE])
    }

    # Stop the parallel backend
    stopImplicitCluster()
    stopCluster(cl)
    # --------------------

    # Clear some memory
    gc()

    if (i_main == 1) {
      if (calculate_silhoutte | calculate_davies_bouldin_index) {
        no_factor_cluster <- as.numeric(levels(initial_clustering))[initial_clustering]  # This weird line converts a vector with factors to a normal numeric vector.
      }
      if (calculate_silhoutte) {
        print(paste("  Calculating silhoutte of likelihood-clusters for iteration 1"), quote = FALSE)
        silhouette_measure <- mean(cluster::silhouette(no_factor_cluster, dist(normalized_posterior_probability))[, 3])
      }else {
        silhouette_measure <- 0
      }
      if (calculate_davies_bouldin_index) {
        print(paste("  Calculating Davies Bouldin's index of likelihood-clusters for iteration 1"), quote = FALSE)
        db <- clusterSim::index.DB(normalized_posterior_probability, no_factor_cluster)$DB
      }else {
        db <- 0
      }
    }


    ####################################################################
    ##### C-step #######################################################
    ##### update cluster allocations ###################################
    print(paste("  Assigning new clusters"), quote = FALSE)

      updated_cell_clust <- sapply(seq_len(nrow(normalized_posterior_probability)),
                                   function(row) which.max(normalized_posterior_probability[row,]))

    updated_cell_clust <- unlist(updated_cell_clust)

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
        converged <- T
        stop_iterating_flag <- T
        break
      }
    }

    rand_index_true_cluster <- aricode::RI(true_cell_cluster_allocation, updated_cell_clust)
    RIs_all[i_main] <- rand_index_true_cluster
    time_taken <- hms_span(start_time_iteration, Sys.time())
    print(paste0(" Iteration ", i_main, ", took ", time_taken, ", Rand index: ", rand_index_true_cluster), quote = FALSE)


    if(calculate_optimization_target){
      # Calculate max loglikeihood from eq 2.2 in the CEM paper
      # rowSums(unnormalized_joint_posterior) is n_cells long, max_loglikelihood is a scalar
      rowsums_unnormalized_joint_posterior <- rowSums(unnormalized_joint_posterior)
      max_loglikelihood <- Rmpfr::asNumeric(sum(log(rowsums_unnormalized_joint_posterior)))  # according to equation 2.2 in CEM paper
      max_loglikelihood_all[i_main] <- max_loglikelihood
      # Calculate CML_C2 from eq. 2.4 in CEM paper
      CML_C2 <- as.numeric(as.character(current_cell_cluster_allocation))
      CML_C2 <- Rmpfr::mpfr(x = CML_C2, precBits = 512)
      for(i_row in seq_along(current_cell_cluster_allocation)){
        c_col <- current_cell_cluster_allocation[i_row]
        CML_C2[i_row] <- log(unnormalized_joint_posterior[i_row, c_col])
      }
      CML_C2_all[i_main] <- Rmpfr::asNumeric(sum(CML_C2))
    }

    current_cell_cluster_allocation <- as.factor(updated_cell_clust)
    if (stop_iterating_flag) {
      # Clean up cluster history
      cell_cluster_history <- cell_cluster_history[, colSums(is.na(cell_cluster_history)) == 0, drop = FALSE]
      # Stop iterations/exit function
      break
    }

    print(paste("  Plotting for iteration."), quote = FALSE)

    # scatter_plot_loglikelihood(dat = dat,
    #                            likelihood = normalized_posterior_probability,
    #                            n_cell_clusters = n_cell_clusters,
    #                            penalization_lambda = penalization_lambda,
    #                            output_path = output_path,
    #                            i_main = i_main,
    #                            true_cell_cluster_allocation_vector = true_cell_cluster_allocation)


    # hist_plot_loglikelihood(dat = dat,
    #                         likelihood = normalized_posterior_probability,
    #                         n_cell_clusters = n_cell_clusters,
    #                         penalization_lambda = penalization_lambda,
    #                         output_path = output_path,
    #                         i_main = i_main,
    #                         true_cell_cluster_allocation_vector = true_cell_cluster_allocation)
  }

  start_time_alluvial_plot <- Sys.time()

  cell_cluster_history_plotting <- cbind('Cell ID' = cell_cluster_history[, 1],
                                         'True cell\ncluster allocation' = true_cell_cluster_allocation,
                                         cell_cluster_history[, 2:ncol(cell_cluster_history)])
  # png(file.path(output_path, paste0("Alluvial_diagram_lambda_", round(penalization_lambda, 6),"_", plot_suffix, ".png")),
  #     width = 1024 + ncol(cell_cluster_history_plotting) * 40, height = 1024, units = "px", res = 150)
  pdf(file.path(output_path, paste0("Alluvial_diagram_lambda_", round(penalization_lambda, 6),"_", plot_suffix, ".pdf")),
      width = 6  + ncol(cell_cluster_history_plotting) * 0.5, height = 6)
  plot_cluster_history(cell_cluster_history = cell_cluster_history_plotting, correct_plot = FALSE)
  dev.off()
  time_taken <- round(Sys.time() - start_time_alluvial_plot, 2)
  print(paste(" Iterations complete, Alluvial plot took", time_taken, "seconds."), quote = FALSE)

  return(
    list(
      "cell_cluster_allocation" = cell_cluster_history[,ncol(cell_cluster_history)][[1]],
      "scregclust_final_result_models" = scregclust_final_result_models,
      "scregclust_final_result_module" = scregclust_final_result_module,
      "scregclust_final_result_coeffs" = scregclust_final_result_coeffs,
      "scregclust_final_result_signs" = scregclust_final_result_signs,
      "scregclust_final_result_reg_table" = scregclust_final_result_reg_table,
      "rand_index" = rand_index_true_cluster,
      "n_iterations" = i_main,
      "silhouette_measure" = silhouette_measure,
      "davies_bouldin_index" = db,
      "RIs_all" = RIs_all[1:i_main],
      "max_loglikelihood" = max_loglikelihood_all[1:i_main],
      "CML_C2" = CML_C2_all[1:i_main],
      "target_genes_residual_var" = target_genes_residual_var,
      "converged" = converged,
      "call" = list(
        "cell_id" = cell_id,
        "true_cell_cluster_allocation" = true_cell_cluster_allocation,
        "max_iter" = max_iter,
        "n_target_gene_clusters" = n_target_gene_clusters,
        "initial_clustering" = initial_clustering,
        "n_cell_clusters" = n_cell_clusters,
        "ind_targetgenes" = ind_targetgenes,
        "ind_reggenes" = ind_reggenes,
        "output_path" = output_path,
        "penalization_lambda" = penalization_lambda,
        "calculate_optimization_target" = calculate_optimization_target,
        "calculate_silhoutte" = calculate_silhoutte,
        "calculate_davies_bouldin_index" = calculate_davies_bouldin_index,
        "plot_suffix" = plot_suffix,
        "always_use_flat_prior" = always_use_flat_prior,
        "use_garbage_cluster_targets" = use_garbage_cluster_targets,
        "retain_gene_clusters" = retain_gene_clusters
      ),
      "metaparameters" = list(
        "normalisation_parameters" = normalisation_parameters,          # parameters necessary to recreate normalization function
        "likelihood_models"        = models                             # neat coefficients because I don't want to sort them out again
      )
    )
  )
}


# Runs only when script is run by itself
# || interactive()
if (sys.nframe() == 0) {

}
