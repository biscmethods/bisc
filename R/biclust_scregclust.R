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


  for (i_cell_cluster in seq_len(n_cell_clusters)) {
    # Standardize data like in scregclust
    cell_cluster_rows_ind <- which(current_cell_cluster_allocation == i_cell_cluster)
    training_data_ind <- which(data_split_for_scregclust[[i_cell_cluster]] == 1)
    training_data_for_scregclust_ind <- cell_cluster_rows_ind[training_data_ind]
    rest_of_data_ind <- which(!inverse_which(training_data_for_scregclust_ind, n_cells))
    res <- standardize_like_scregclust(xvals = regulator_genes, yvals = target_genes,
                                       training_data_ind = training_data_for_scregclust_ind, test_data_ind = rest_of_data_ind)
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
  return(loglikelihood)
}


standardize_like_scregclust <- function(xvals, yvals, training_data_ind, test_data_ind) {


  org_yvals <- yvals
  org_xvals <- xvals

  yvals[training_data_ind,] <- scale(yvals[training_data_ind,], scale = FALSE)
  xvals[training_data_ind,] <- scale(xvals[training_data_ind,])

  yvals[test_data_ind,] <- scale(yvals[test_data_ind,], colMeans(org_yvals[training_data_ind,]), scale = FALSE)
  xvals[test_data_ind,] <- scale(xvals[test_data_ind,], colMeans(org_xvals[training_data_ind,]), apply(org_xvals[training_data_ind,], 2, sd))

  return(list("xvals" = xvals, "yvals" = yvals))
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

        if ((length(ind_reggenes) > n_test_data) || (length(ind_reggenes) > n_training_data)) {
          print(paste("   The split into training/test data for scregclust resulted in less cells than there are regulators."), quote = FALSE)
          print(paste("    There are", length(ind_reggenes), "regulator genes."), quote = FALSE)
          print(paste0("    There are ", ncol(indata_for_scregclust), " cells in cell cluster ", i_cell_cluster, "."), quote = FALSE)
          print(paste("    There are", n_test_data, "cells in the test data."), quote = FALSE)
          print(paste("    There are", n_training_data, "cells in the training data."), quote = FALSE)
          return(NA)
        }
        sink(tempfile()) # Shut scregclust up
        screg_out <- scregclust::scregclust(
          expression = indata_for_scregclust,  # scRegClust wants this form
          split_indices = data_split_for_scregclust[[i_cell_cluster]],
          genesymbols = 1:(n_target_genes + n_regulator_genes),  # Gene row numbers
          is_regulator = inverse_which(indices = ind_reggenes, output_length = n_regulator_genes + n_target_genes),  # Vector indicating which genes are regulators
          n_cl = n_target_gene_clusters[[i_cell_cluster]],
          penalization = penalization_lambda,
          noise_threshold = 0.00001,
          verbose = FALSE,
          # n_cycles = 100,
        )
        sink()

        # This is OK, things work out anyway
        # if (!screg_out$results[[1]]$converged) {
        #   print(paste("scregclust didn't converge for cell cluster", i_cell_cluster,"."), quote = FALSE)
        #   return(NA)
        # }

        screg_out_betas <- do.call(cbind, screg_out$results[[1]]$output[[1]]$coeffs)  # Merge betas into one matrix
        if (is.null(screg_out_betas)) {
          print(paste("Cell cluster", i_cell_cluster, "betas is NULL for scregclust output."), quote = FALSE)
          return(NA)
        }

        # Create cluster_indexes which maps columns in screg_out_betas to correct places under some assumptions of order:
        # target_gene_cluster_names, e.g. 2 2 1 3 1 3,
        # becomes cluster_indexes: 3 4 1 5 2 6
        # CAN ALSO BE DONE WITH THE FOLLOWING CODE THAT WE HAVE USED BEFORE:
        # indexes_of_not_deleted_target_gene_clusters <- which(sapply(out_list[[i_target_cell_cluster]]$results[[1]]$output[[1]]$coeffs, FUN=function(x) !is.null(x)))
        # target_gene_cluster_index <- out_list[[i_target_cell_cluster]]$results[[1]]$output[[1]]$cluster[1:n_target_genes]  # n_target_genes, e.g. 1 1 1 1 3 3 3 3 3 3 3 3 1 3 1 3 3 3 3 3 1 3 1 1 3 1 3 3 1 3
        # target_gene_cluster_index <- unlist(sapply(indexes_of_not_deleted_target_gene_clusters, FUN=function(x) which(target_gene_cluster_index==x)))
        target_gene_cluster_names <- screg_out$results[[1]]$output[[1]]$cluster[1:n_target_genes]
        cluster_indexes <- rep(NA, length(target_gene_cluster_names))
        tmp_max_ind <- 0
        for (i_target_gene_cluster in 1:n_target_gene_clusters[[i_cell_cluster]]) {
          ind_tmp <- which(target_gene_cluster_names == i_target_gene_cluster)
          n_target_genes_in_cluster <- length(ind_tmp)
          cluster_indexes[ind_tmp] <- 1:n_target_genes_in_cluster + tmp_max_ind  # Get names of clusters for each target gene
          tmp_max_ind <- n_target_genes_in_cluster
        }

        screg_out_betas <- screg_out_betas[, cluster_indexes]  # Sort the matrix


      }
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

    cluster_proportions <- vector(length = n_cell_clusters)
    for (i_cell_cluster in seq_len(n_cell_clusters)) {
      print(paste("  Calculating cluster proportions for cell cluster", i_cell_cluster), quote = FALSE)
      current_cluster_proportion <- sum(current_cell_cluster_allocation == i_cell_cluster) / length(current_cell_cluster_allocation)

      # If some cluster has been completely omitted, give it a nonzero proportion anyway for next iteration
      # Even without weights this is good to include for cluster allocation (though it is a bit arbitrary)
      if (current_cluster_proportion == 0) {
        cluster_proportions[i_cell_cluster] <- 0.0001
        stop_iterating_flag <- T
      }else
        cluster_proportions[i_cell_cluster] <- current_cluster_proportion
    }

    ####################################################################
    ##### E-step #######################################################
    ####### Compute posterior probabilities ############################
    ####### For each cell to belong to each cluster ####################

    # If everything goes well this should work the same for screg as with our
    # lm examples
    # weights <- sweep(exp(loglikelihood/100), 2, cluster_proportions, "*")

    print(paste("  Doing final weights calculations"), quote = FALSE)
    loglikelihood <- Rmpfr::mpfr(x = loglikelihood, precBits = 256)
    weights <- sweep(exp(loglikelihood), 2, (cluster_proportions), "*")

    big_logLhat_in_BIC <- sum(log(rowSums(weights)))
    k_in_BIC <- (n_target_genes + n_regulator_genes) * n_cell_clusters
    n_in_BIC <- n_total_cells
    BIC <- k_in_BIC * log(n_in_BIC) - 2 * big_logLhat_in_BIC

    weights <- sweep(weights, 1, rowSums(weights), "/")
    weights <- Rmpfr::asNumeric(weights)  # if you use nórmal as.numeric it doens't keep the matrix format

    loglikelihood <- weights

    weights_all[[i_main]] <- weights
    BIC_all[[i_main]] <- BIC

    likelihood_all[[i_main]] <- loglikelihood


    if (i_main == 1) {
      print(paste("  Calculating complexity of likelihood-clusters for iteration 1"), quote = FALSE)
      no_factor_cluster <- as.numeric(levels(initial_clustering))[initial_clustering]  # This weird line converts a vector with factors to a normal numeric vector.
      # print(str(no_factor_cluster), quote = FALSE)
      # db <- index.DB(loglikelihood, no_factor_cluster)$DB
      db <- mean(silhouette(no_factor_cluster, dist(loglikelihood))[, 3])
    }


    ####################################################################
    ##### C-step #######################################################
    ##### update cluster allocations ###################################
    print(paste("  Assigning new clusters"), quote = FALSE)
    if (use_complex_cluster_allocation) {
      pca_likelihood <- stats::prcomp(loglikelihood, center = TRUE, scale. = TRUE)$x
      fcm_res <- ppclust::fcm(x = pca_likelihood, centers = n_cell_clusters)
      weights_all[[i_main]] <- fcm_res$u
      updated_cell_clust <- fcm_res$cluster
      # updated_cell_clust <- stats::kmeans(x = pca_likelihood, centers = n_cell_clusters, iter.max = 20, nstart = 50 + n_cell_clusters)$cluster
    }else {
      updated_cell_clust <- sapply(seq_len(nrow(loglikelihood)),
                                   function(row) which.max(loglikelihood[row,]))
    }
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
        stop_iterating_flag <- T
        break
      }
    }

    rand_index_true_cluster <- aricode::RI(true_cell_cluster_allocation, updated_cell_clust)
    time_taken <- hms_span(start_time_iteration, Sys.time())
    print(paste0(" Iteration ", i_main, ", took ", time_taken, ", Rand index: ", rand_index_true_cluster), quote = FALSE)

    current_cell_cluster_allocation <- as.factor(updated_cell_clust)

    if (stop_iterating_flag) {
      # Clean up cluster history
      cell_cluster_history <- cell_cluster_history[, colSums(is.na(cell_cluster_history)) == 0, drop = FALSE]
      # Stop iterations/exit function
      break
    }

    print(paste("  Plotting for iteration."), quote = FALSE)

    scatter_plot_loglikelihood(dat = dat,
                               likelihood = loglikelihood,
                               n_cell_clusters = n_cell_clusters,
                               penalization_lambda = penalization_lambda,
                               output_path = output_path,
                               i_main = i_main,
                               true_cell_cluster_allocation_vector = true_cell_cluster_allocation)

    hist_plot_loglikelihood(dat = dat,
                            likelihood = loglikelihood,
                            n_cell_clusters = n_cell_clusters,
                            penalization_lambda = penalization_lambda,
                            output_path = output_path,
                            i_main = i_main,
                            true_cell_cluster_allocation_vector = true_cell_cluster_allocation)
  }

  start_time_alluvial_plot <- Sys.time()

  cell_cluster_history_plotting <- cbind('Cell ID' = cell_cluster_history[, 1],
                                         'True cell cluster allocation' = true_cell_cluster_allocation,
                                         cell_cluster_history[, 2:ncol(cell_cluster_history)])
  png(file.path(output_path, paste0("Alluvial_diagram_lambda_", round(penalization_lambda, 6), ".png")),
      width = 1024 + ncol(cell_cluster_history_plotting) * 40, height = 1024, units = "px", res = 150)
  plot_cluster_history(cell_cluster_history = cell_cluster_history_plotting, correct_plot = FALSE)
  dev.off()
  time_taken <- round(Sys.time() - start_time_alluvial_plot, 2)
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

}
