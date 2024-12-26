# !/usr/bin/Rscript
# # # for dev
#  new_data <- scenarios[[1]]$biclust_input_data[1:10,]
# #
#  fitted_model <- bisc_results_list[[1]]$bisc_results[[1]]
#
#  bisc_results_list[[1]]$bisc_results[[1]]$target_genes_residual_var
#
# #
#
# prior_cluster_proportions = NULL
# calculate_optimization_target             = TRUE
# seed = 1234
#
# #
# str(bisc_results_list[[1]]$bisc_results[[1]]$scregclust_final_result_coeffs)
# bisc_results_list[[1]]$bisc_results[[1]]$call$ind_reggenes
# bisc_results_list[[1]]$bisc_results[[1]]$metaparameters$normalisation_parameters


# function to predict cell cluster allocations of new data given a previously fitted model
# returns new cell cluster allocation

library(here)
R_path <- here::here("R")
# source(file.path(R_path, "bisc.R"))


#local standardisation
standardize_like_scregclust_ <- function(xvals, yvals,
                                         normalisation_parameters
) {

  yval_colmeans  <- normalisation_parameters$yval_colmeans
  xvals_colmeans <- normalisation_parameters$xvals_colmeans
  xval_scale     <- normalisation_parameters$xval_scale

  yvals[,] <- scale(yvals[,], yval_colmeans, scale = FALSE)
  xvals[,] <- scale(xvals[,], xvals_colmeans, xval_scale )

  return(list("xvals" = xvals, "yvals" = yvals,
              "params" = list(
                'yval_colmeans'  = yval_colmeans,
                'xvals_colmeans' = xvals_colmeans,
                'xval_scale'     = xval_scale
              )
  )
  )
}

# local likelihood calc

loglikelihood_calc_matrix_ <- function(all_genes, # Data
                                       models,
                                       target_genes_residual_var,
                                       n_cell_clusters,
                                       ind_reggenes,
                                       ind_targetgenes,
                                       # data_split_for_scregclust,
                                       # current_cell_cluster_allocation,
                                       normalisation_parameters
) {
  # For all cells, calculate the likelihood of coming from the model corresponding to each
  n_cells <- nrow(all_genes)
  loglikelihood <- matrix(data = 0, nrow = n_cells, ncol = n_cell_clusters)

  regulator_genes <- as.matrix(all_genes[, ind_reggenes, drop = FALSE])  # 1xr -> cxr
  target_genes <- as.matrix(all_genes[, ind_targetgenes, drop = FALSE])  # 1xt -> cxt

  parameters <- vector(mode = "list", length = n_cell_clusters)

  for (i_cell_cluster in seq_len(n_cell_clusters)) {
    # Standardize data like in scregclust
    # cell_cluster_rows_ind <- which(current_cell_cluster_allocation == i_cell_cluster)

    # training_data_ind <- which(data_split_for_scregclust[[i_cell_cluster]] == 1)

    # training_data_for_scregclust_ind <- cell_cluster_rows_ind[training_data_ind]

    # rest_of_data_ind <- which(!inverse_which(training_data_for_scregclust_ind, n_cells))

    #save this function and return it for use in the predict function

    res <- standardize_like_scregclust_(xvals = regulator_genes, yvals = target_genes,
                                        normalisation_parameters = normalisation_parameters[[i_cell_cluster]]
    )

    parameters[[i_cell_cluster]] <- res$params

    regulator_genes_normalized <- res[['xvals']]

    target_genes_normalized <- res[['yvals']]

    #--------------------

    print(paste("  Calculating loglikelihood for cell cluster", i_cell_cluster), quote = FALSE)


    cell_cluster_betas <- models[[i_cell_cluster]]     # rxt

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
    'parameters' = res$params
  )
  )
}






bisc_predict <- function(new_data,     # as in scenarios[[1]]$biclust_input_data
                         fitted_model, # as in bisc_results_list[[1]]$bisc_results[[1]]
                         prior_cluster_proportions = NULL,
                         calculate_optimization_target             = TRUE,

                         seed = 1234){

  set.seed(seed)

  # extract parameters from fit
  models                          <- fitted_model$metaparameters$likelihood_models
  target_genes_residual_var       <- fitted_model$target_genes_residual_var
  n_cell_clusters                 <- fitted_model$call$n_cell_clusters
  ind_reggenes                    <- fitted_model$call$ind_reggenes
  ind_targetgenes                 <- fitted_model$call$ind_targetgenes
  current_cell_cluster_allocation <- fitted_model$cell_cluster_allocation
  n_total_cells                   <- nrow(new_data)
  n_target_genes                  <- length(ind_targetgenes)
  n_regulator_genes               <- length(ind_reggenes)

  # calculate necessary loglikelihoods for new data
  loglikelihood <- loglikelihood_calc_matrix_(all_genes = new_data,
                                              models = models,
                                              target_genes_residual_var = target_genes_residual_var,
                                              n_cell_clusters = n_cell_clusters,
                                              ind_reggenes = ind_reggenes,
                                              ind_targetgenes = ind_targetgenes,
                                              normalisation_parameters = fitted_model$metaparameters$normalisation_parameters
  )

  # separate out normalisation parameters
  normalisation_parameters <- loglikelihood$parameters
  loglikelihood            <- loglikelihood$loglikelihood


  if (!is.null(prior_cluster_proportions)) {

    cluster_proportions <- prior_cluster_proportions

    # cluster_proportions <- vector(length = n_cell_clusters)
    #
    # for (i_cell_cluster in seq_len(n_cell_clusters)) {
    #   print(paste("  Calculating cluster proportions for cell cluster", i_cell_cluster), quote = FALSE)
    #   current_cluster_proportion <- sum(current_cell_cluster_allocation == i_cell_cluster) / length(current_cell_cluster_allocation)
    #
    #   # If some cluster has been completely omitted, give it a nonzero proportion anyway for next iteration
    #   # Even without weights this is good to include for cluster allocation (though it is arbitrary)
    #   if (current_cluster_proportion == 0) {
    #     cluster_proportions[i_cell_cluster] <- 0.0001
    #     stop_iterating_flag <- T
    #   }else
    #     cluster_proportions[i_cell_cluster] <- current_cluster_proportion
    # }

  } else {
    cluster_proportions <- (vector(length = n_cell_clusters) + 1) / n_cell_clusters
  }


  ####################################################################
  ##### E-step #######################################################
  ####### Compute posterior probabilities ############################
  ####### For each cell to belong to each cluster ####################

  # If everything goes well this should work the same for screg as with our
  # lm examples


  print(paste("  Doing final weights calculations"), quote = FALSE)
  loglikelihood <- Rmpfr::mpfr(x = loglikelihood, precBits = 256)

  weights <- sweep(exp(loglikelihood), 2, (cluster_proportions), "*")


  if (calculate_optimization_target) {
    big_logLhat_in_BIC <- sum(log(rowSums(weights)))
    k_in_BIC <- (n_target_genes + n_regulator_genes) * n_cell_clusters
    n_in_BIC <- n_total_cells
    BIC <- k_in_BIC * log(n_in_BIC) - 2 * big_logLhat_in_BIC
  }else {
    BIC <- 0
  }

  # Below parallizes the this calculation (e.g speeds up from 1:32 min to 12 sec)
  # weights <- sweep(weights, 1, rowSums(weights), "/")

  # Define the number of cores to use
  num_cores <- max(detectCores() - 2, 1)

  # # Create a parallel backend
  cl <- parallel::makeCluster(num_cores,  type = "PSOCK")

  # Register the parallel backend
  registerDoParallel(cl)
  num_rows <- nrow(weights)

  split_indices <- split(1:num_rows, cut(1:num_rows, num_cores))

  # Define the function to parallelize
  parallel_sweep <- function(x) {
    Rmpfr::asNumeric(sweep(x, 1, rowSums(x), "/"))
  }

  # Parallelize the operation using foreach
  weights <- foreach(i = split_indices, .packages = c("Rmpfr"), .combine = 'rbind') %dopar% {
    parallel_sweep(weights[i, , drop = FALSE])
  }

  # Stop the parallel backend
  stopImplicitCluster()
  stopCluster(cl)


  # --------------------

  loglikelihood <- weights

  ####################################################################
  ##### C-step #######################################################
  ##### update cluster allocations ###################################

  #requires: loglikelihoods
  print(paste("  Assigning new clusters"), quote = FALSE)
  updated_cell_clust <- sapply(seq_len(nrow(loglikelihood)),
                               function(row) which.max(loglikelihood[row,]))
  updated_cell_clust <- unlist(updated_cell_clust)



  return(
    list(
      "cell_cluster_allocation" = updated_cell_clust
    )
  )
}


# #dev test
#
# out <- bisc_predict(new_data,     # as in scenarios[[1]]$biclust_input_data
#              fitted_model, # as in bisc_results_list[[1]]$bisc_results[[1]]
#              prior_cluster_proportions = NULL,
#              calculate_optimization_target             = TRUE,
#
#              seed = 1234)

