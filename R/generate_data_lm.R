library(tidyverse)

#' Helper function for checking arguments
#'
#' Checks if an input argument is e.g. positive.
#' @param x The input argument you want to check.
#' @param one_element Checks if it's just one element.
#' @param atomic Checks if it's atomic, i.e. it's not a list.
#' @param numeric Checks if it's numeric.
#' @param positive Checks if it's positive.
#' @param int Checks if it's an integer.
#' @return void, or throws an error message
#' @examples
#'   checks(x=1, one_element=TRUE, atomic=TRUE, numeric=TRUE, positive=TRUE, int=FALSE)
#' @noRd
checks <- function(x,
                   one_element = TRUE,
                   atomic = TRUE,
                   numeric = TRUE,
                   positive = TRUE,
                   int = TRUE) {
  error_message <- ""
  if (one_element) {
    if (!(length(x) == 1L)) {
      add_message <- paste0(deparse(substitute(x)), " must be one number.")
      error_message <- paste(error_message, add_message, sep = "\n")
    }
  }

  if (atomic) {
    if (!is.atomic(x)) {
      add_message <- (paste0(deparse(substitute(x)), " must be atomic."))
      error_message <- paste(error_message, add_message, sep = "\n")
    }
  }

  if (numeric) {
    if (!is.numeric(x)) {
      add_message <- (paste0(deparse(substitute(x)), " must be numeric."))
      error_message <- paste(error_message, add_message, sep = "\n")
    }
  }

  if (positive) {
    if (!all(x > 0)) {
      add_message <- (paste0(deparse(substitute(x)), " must be >0."))
      error_message <- paste(error_message, add_message, sep = "\n")
    }
  }

  if (int) {
    if (!(all(round(x) == x))) {
      add_message <- (paste0(deparse(substitute(x)), " must be an integer/s, now it's: ", toString(x)))
      error_message <- paste(error_message, add_message, sep = "\n")
    }
  }

  if (nchar(error_message) > 1) {
    stop(error_message)
  }
}


#' Dummy data generation for biclust with lm
#'
#' Generates dummy data that work with biclust with lm
#' For now we don't generate sparse data with regulators and targets set to 0. Also known as, we don't generate target gene type clusters.
#' That's for later when we want to play with e.g. lasso.
#' We assume:
#' ERROR: The error variance is the same across all cell clusters. The formal assumption for (many x -> many y linear regression) is:
#'   The errors are usually assumed to be uncorrelated across measurements, and follow a multivariate normal distribution. If the errors do not follow a multivariate normal distribution, generalized linear models may be used to relax assumptions about Y and U.
#' BETAS/COEFFICIENTS' MEANS AND SD:
#'   Close to zero mean and high variance, one different mean and sd per target gene (same across all regulator genes). This arranges it so that we have one slope per cell cluster.
#' REGULATOR MEANS AND SD:
#'   One mean and sd per cell cluster. It might be easier just to keep it the same for all cell clusters.
#' TODO: Go through the theory and assumptions of the model: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
#' @param n_cell_clusters The number of cell clusters.
#' @param n_target_gene_type The number of target gene types. We have x named target genes that have one expression per cell.
#' @param n_regulator_gene_type  The number of regulator gene types. We have x named regulator genes that have one expression per cell.
#' @param n_cells c(1000,5000,10000), The number of cells in each cell cluster given as a vector.
#' @param regulator_means c(1,2,3), Regulator gene expression mean in each cell cluster.
#' @param regulator_standard_deviations = c(0.1,0.2,0.3),  Regulator sd for expression in each cell cluster.
#' @param coefficients_standard_deviation = 100, 'The betas/slopes'. One per target gene. Instead of providing mean and standard deviation for each target gene, we provide the standard deviation from which these will be generated. Mean will be 0.
#' @param target_gene_type_standard_deviation = 3, The error, the sd of the expression of the target genes.
#' @return dat
#' @examples
#'   dat <- generate_data_lm(n_cell_clusters = 3,
#'                           n_target_gene_type = 5,
#'                           n_regulator_gene_type = 20,
#'                           n_cells = c(1000,5000,10000),
#'                           regulator_means = c(1,2,3),
#'                           regulator_standard_deviations = c(0.1,0.2,0.3),
#'                           coefficients_standard_deviation = 100,
#'                           target_gene_type_standard_deviation = 3,
#'                           )
#' @export
generate_data_lm <- function(n_cell_clusters = 3,
                             n_target_gene_type = 5,  # We have x named target genes that have one expression per cell
                             n_regulator_gene_type = 20,  # We have x named regulator genes that have one expression per cell
                             n_cells = c(1000, 5000, 10000),
                             regulator_means = c(1, 2, 3),  # Regulator mean expression in each cell cluster.
                             regulator_standard_deviations = c(0.1, 0.2, 0.3),  # Regulator sd for expression in each cell cluster.
                             coefficients_standard_deviation = 100, # 'The betas/slopes'. One per target gene. Instead of providing mean and standard deviation for each target gene, we provide the standard deviation from which these will be generated. Mean will be 0.
                             target_gene_type_standard_deviation = 3
) {


  # Manual arguments ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # n_cell_clusters <- 3
  # n_target_gene_type <- 5  # We have x named target genes that have one expression per cell
  # n_regulator_gene_type <- 20  # We have x named regulator genes that have one expression per cell
  # n_cells <- c(1000, 5000, 10000)
  # regulator_means <- c(1, 2, 3)  # For generating dummy data, regulator mean in each cell cluster
  # regulator_standard_deviations <- c(0.1, 0.2, 0.3)  # For generating dummy data, regulator mean in each cell cluster
  # coefficients_standard_deviation <- 100 # 'The betas/slopes'. One per target gene. Instead of providing mean and standard deviation for each target gene, we provide the standard deviation from which these will be generated. Mean will be 0.
  # target_gene_type_standard_deviation <- 3


  # Check arguments -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  checks(n_cell_clusters, one_element = TRUE, atomic = TRUE, numeric = TRUE, positive = TRUE, int = TRUE)
  checks(n_target_gene_type, one_element = TRUE, atomic = TRUE, numeric = TRUE, positive = TRUE, int = TRUE)
  checks(n_regulator_gene_type, one_element = TRUE, atomic = TRUE, numeric = TRUE, positive = TRUE, int = TRUE)
  checks(n_cells, one_element = FALSE, atomic = TRUE, numeric = TRUE, positive = TRUE, int = TRUE)
  checks(regulator_means, one_element = FALSE, atomic = TRUE, numeric = TRUE, positive = FALSE, int = FALSE)
  checks(regulator_standard_deviations, one_element = FALSE, atomic = TRUE, numeric = TRUE, positive = TRUE, int = FALSE)
  checks(coefficients_standard_deviation, one_element = TRUE, atomic = TRUE, numeric = TRUE, positive = TRUE, int = FALSE)
  checks(target_gene_type_standard_deviation, one_element = TRUE, atomic = TRUE, numeric = TRUE, positive = TRUE, int = FALSE)

  if (length(n_cells) != n_cell_clusters) {
    stop(paste0("n_cells has length ", length(n_cells), ", but must have length n_cell_clusters: ", n_cell_clusters, "."))
  }

  if (length(regulator_means) != n_cell_clusters) {
    stop(paste0("regulator_means has length ", length(regulator_means), ", but must have length n_cell_clusters: ", n_cell_clusters, "."))
  }

  if (length(regulator_standard_deviations) != n_cell_clusters) {
    stop(paste0("regulator_standard_deviations has length ", length(regulator_standard_deviations), ", but must have length n_cell_clusters: ", n_cell_clusters, "."))
  }


  # Calculated variables --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  n_total_cells <- sum(n_cells)
  true_cell_cluster_allocation <- rep(1:n_cell_clusters, times = n_cells)  # Cell cluster IDs will be in order


  # Generate regulator expressions for each cell cluster ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  regulator_gene_expression <- vector(mode = 'list', length = n_cell_clusters)
  for (i_cell_cluster in 1:n_cell_clusters) {
    current_n_cells <- n_cells[i_cell_cluster]
    temp_regulator_gene_expressions_vector <- rnorm(current_n_cells * n_regulator_gene_type,
                                                    mean = regulator_means[i_cell_cluster],
                                                    sd = regulator_standard_deviations[i_cell_cluster])
    regulator_gene_expression[[i_cell_cluster]] <- matrix(temp_regulator_gene_expressions_vector,
                                                          ncol = n_regulator_gene_type)
  }
  regulator_gene_expression <- do.call(rbind, regulator_gene_expression)
  colnames(regulator_gene_expression) <- paste0("r", 1:n_regulator_gene_type)
  regulator_gene_expression <- tibble::as_tibble(regulator_gene_expression)


  # Generate betas --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # One per target gene type in each cell cluster meaning:
  # Two cell clusters with one target gene type in each: total number of betas=2
  # Two cell clusters with 3 respectively 4 target gene types: total number of betas=7
  betas <- vector(mode = 'list', length = n_cell_clusters)
  for (i_cell_cluster in 1:n_cell_clusters) {
    temp_betas_for_cell_cluster <- rnorm(n_target_gene_type,
                                         mean = 0,
                                         sd = coefficients_standard_deviation)
    # We want the different betas for each target gene type in the cell cluster, but the same across all regulator genes
    # Produces a matrix (Target gene types) x (Regulator gene types)
    temp_betas <- t(matrix(rep(temp_betas_for_cell_cluster, n_regulator_gene_type), ncol = n_regulator_gene_type))
    colnames(temp_betas) <- paste0("b", 1:n_target_gene_type)
    betas[[i_cell_cluster]] <- tibble::as_tibble(temp_betas)
  }


  # Calculate target gene expressions for each cell cluster ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  cell_cluster_target_gene_expression <- sapply(1:n_cell_clusters,
                                                function(i_cell_cluster) {
                                                  current_betas <- betas[[i_cell_cluster]]  # This picks out the relevant betas, as many as there are target gene types in this cell cluster
                                                  current_regulator_gene_index <- which(true_cell_cluster_allocation == i_cell_cluster)
                                                  current_regulator_gene_expressions <- regulator_gene_expression[current_regulator_gene_index,]
                                                  residuals <- rnorm(n_cells[i_cell_cluster],
                                                                     mean = 0,
                                                                     sd = target_gene_type_standard_deviation)
                                                  return(as.matrix(current_regulator_gene_expressions) %*% as.matrix(current_betas) + residuals)
                                                }
  )
  # Put it into a tibble
  target_gene_expression <- do.call(rbind, cell_cluster_target_gene_expression)
  colnames(target_gene_expression) <- paste0("t", 1:n_target_gene_type)
  target_gene_expression <- tibble::as_tibble(target_gene_expression)


  # Construct output data structure (also a tibble) -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  dat <- tibble::tibble(cell_id = (1:n_total_cells),
                        true_cell_cluster_allocation = factor(true_cell_cluster_allocation),
                        target_gene_expression,
                        regulator_gene_expression)

  return(dat)
  # return indexes of different gene types?
  # list(dat = dat,
  #      Z_r = Z_r,
  #      Pi  = Pi,
  #      R   = R,
  #      S   = S,
  #      B = Beta_with_signs
  # )
}


# Runs only when script is run by itself
# || interactive()
if (sys.nframe() == 0) {
  # ... do main stuff
  library(ggfortify)  # For pca-plot

  # Set seed for example
  set.seed(1234)

  dat <- generate_data_lm(n_cell_clusters = 3,
                          n_target_gene_type = 5,  # We have x named target genes that have one expression per cell
                          n_regulator_gene_type = 20,  # We have x named regulator genes that have one expression per cell
                          n_cells = c(1000, 5000, 10000),
                          regulator_means = c(1, 2, 5),  # Regulator mean expression in each cell cluster.
                          regulator_standard_deviations = c(0.1, 0.2, 0.3),  # Regulator sd for expression in each cell cluster.
                          coefficients_standard_deviation = 100, # 'The betas/slopes'. One per target gene. Instead of providing mean and standard deviation for each target gene, we provide the standard deviation from which these will be generated. Mean will be 0.
                          target_gene_type_standard_deviation = 3
  )

  dat[, 'true_cell_cluster_allocation'] <- paste("Cluster", pull(dat, 'true_cell_cluster_allocation'))  # These needs to be strings for discrete labels in pca plot
  pca_res <- prcomp(dat[, 3:ncol(dat)], scale. = TRUE)
  p <- ggplot2::autoplot(pca_res, data = dat, colour = 'true_cell_cluster_allocation')
  plot(p)
  print(dat)
}
