# !/usr/bin/Rscript
library(tidyverse)
library(ggmulti)
library(ggplot2)
library(dplyr)
library(ggfortify)

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
#' @param n_cell_clusters The number of cell clusters.
#' @param n_target_gene_type The number of target gene types. We have x named target genes that have one expression per cell.
#' @param n_regulator_gene_type  The number of regulator gene types. We have x named regulator genes that have one expression per cell.
#' @param n_cells c(1000,5000,10000), The number of cells in each cell cluster given as a vector.
#' @param regulator_means c(1,2,3), Regulator gene expression mean in each cell cluster.
#' @param regulator_standard_deviations = c(0.1,0.2,0.3),  Regulator sd for expression in each cell cluster.
#' @param coefficients_standard_deviation = 100, 'The betas/slopes'. One per target gene. Instead of providing mean and standard deviation for each target gene, we provide the standard deviation from which these will be generated. Mean will be 0.
#' @param target_gene_type_standard_deviation = 3, The error, the sd of the expression of the target genes.
#' @return dat
#' @export
generate_data_lm <- function(n_cell_clusters = 3,
                             n_target_gene_type = 5,  # We have x named target genes that have one expression per cell
                             n_regulator_gene_type = 20,  # We have x named regulator genes that have one expression per cell
                             n_cells = c(1000, 5000, 10000),
                             regulator_means = c(1, 2, 3),  # Regulator mean expression in each cell cluster.
                             regulator_standard_deviations = c(0.1, 0.2, 0.3),  # Regulator sd for expression in each cell cluster.
                             coefficients_standard_deviation = 100, # 'The betas/slopes'. One per target gene. Instead of providing mean and standard deviation for each target gene, we provide the standard deviation from which these will be generated. Mean will be 0.
                             target_gene_type_standard_deviation = 3,
                             plot_stuff = F
) {


  # Manual arguments for dev -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # n_cell_clusters <- 3
  # n_target_gene_type <- 5  # We have x named target genes that have one expression per cell
  # n_regulator_gene_type <- 20  # We have x named regulator genes that have one expression per cell
  # n_cells <- c(1000, 5000, 10000)
  # regulator_means <- c(1, 2, 3)  # For generating dummy data, regulator mean in each cell cluster
  # regulator_standard_deviations <- c(0.1, 0.2, 0.3)  # For generating dummy data, regulator mean in each cell cluster
  # coefficients_standard_deviation <- 100 # 'The betas/slopes'. One per target gene. Instead of providing mean and standard deviation for each target gene, we provide the standard deviation from which these will be generated. Mean will be 0.
  # target_gene_type_standard_deviation <- 3
  # Plot_stuff = F


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
  # One column in beta for each target gene. one row for each beta
  # each column contains the linear model for that target gene.
  # To make the data more complex,
  # TODO: duplicate some target genees, so they have the same model and
  # TODO: make it so that different regulators can be shut off in different clusters
  betas <- vector(mode = 'list', length = n_cell_clusters)
  for (i_cell_cluster in 1:n_cell_clusters) {

    #assuming different targets can have different regulatory models this is more accurate
    temp_betas_for_cell_cluster <- rnorm(n_regulator_gene_type * n_target_gene_type,
                                         mean = 0,
                                         sd = coefficients_standard_deviation)

    # We want the different betas for each target gene type in the cell cluster, but the same across all regulator genes
    # Produces a matrix (Target gene types) x (Regulator gene types)

    temp_betas <- matrix(data = temp_betas_for_cell_cluster,
                         nrow = n_regulator_gene_type,
                         ncol = n_target_gene_type)

    colnames(temp_betas) <- paste0("t", 1:n_target_gene_type)
    rownames(temp_betas) <- paste0("r", 1:n_regulator_gene_type)

    betas[[i_cell_cluster]] <- tibble::as_tibble(temp_betas)
  }

  # Calculate target gene expressions for each cell cluster ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # First calculate them normally and then add noise separately
    # For easier future plotting

  predicted_cell_cluster_target_gene_expression <- lapply(1:n_cell_clusters,
                                                          function(i_cell_cluster) {
                                                            current_betas <- betas[[i_cell_cluster]]
                                                            # This picks out the relevant betas, as many as there are target gene types in this cell cluster
                                                            current_regulator_gene_index <- which(true_cell_cluster_allocation == i_cell_cluster)
                                                            current_regulator_gene_expressions <- regulator_gene_expression[current_regulator_gene_index,]
                                                            return(as.matrix(current_regulator_gene_expressions) %*% as.matrix(current_betas))
                                                          }
  )

  # Todo: remove redundant matrix multiplication
  cell_cluster_target_gene_expression <- lapply(1:n_cell_clusters,
                                                function(i_cell_cluster) {
                                                  current_betas <- betas[[i_cell_cluster]]
                                                  # This picks out the relevant betas, as many as there are target gene types in this cell cluster
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

  #also for predicted values
  predicted_target_gene_expression <- do.call(rbind, predicted_cell_cluster_target_gene_expression)
  colnames(predicted_target_gene_expression) <- paste0("t_hat", 1:n_target_gene_type)
  predicted_target_gene_expression <- tibble::as_tibble(predicted_target_gene_expression)

  # Construct output data structure (also a tibble) -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  dat <- tibble::tibble(cell_id = (1:n_total_cells),
                        true_cell_cluster_allocation = factor(true_cell_cluster_allocation),
                        target_gene_expression,
                        regulator_gene_expression)

  if(plot_stuff){
    # Generate plots for visualisation

    demo_path <- here::here("demo")
    R_path <- here::here("R")
    output_path <- demo_path

    vis_dat <- tibble::tibble(cell_id = (1:n_total_cells),
                              true_cell_cluster_allocation = factor(true_cell_cluster_allocation),
                              targets     = target_gene_expression,
                              regulators  = regulator_gene_expression,
                              predictions = predicted_target_gene_expression,
                              residuals   = (target_gene_expression - predicted_target_gene_expression)
    )

    # Visualise the regulator genes, color cells by cluster
    # This could be done using some kind of principal components, but
    #   We generate target genes independently of each other so why bother

    plot_pca <- function(indata     = vis_dat$regulators,
                         cluster_allocations = vis_dat$true_cell_cluster_allocation,
                         plot_title = "Principal components of regulator genes"
    ){



      # Apply PCA
      pca_result <- prcomp(indata, scale. = TRUE)

      # Extract principal components
      pc_df <- as_tibble(pca_result$x[,1:2])

      # Add labels (if you have any)
      pc_df$label <- cluster_allocations

      # Plot the first two principal components
      plot_name <- "regulators_pca.png"

      # png(file.path(demo_path,plot_name))
      pc_df[sample(nrow(pc_df), 1000),] %>%
        ggplot( aes(x = PC1, y = PC2, colour = label)) +
        geom_point(alpha = 0.1) +
        theme_minimal() +
        labs(title = plot_title) +
        scale_color_discrete(name = "Cell cluster")
      # dev.off()
    }

    # Visualise the regulator genes
    png(file.path(output_path, paste0("lm_data_gen_regulators.png")),
        width = 1024, height = 480, units = "px")

    plot_pca(indata     = vis_dat$regulators,
             plot_title = "Principal components of regulator genes") -> p

    print(p)

    dev.off()

    # Visualise the target gene expressions

    png(file.path(output_path, paste0("lm_data_gen_targets.png")),
        width = 1024, height = 480, units = "px")

    plot_pca(indata = vis_dat$targets,
             plot_title = "Principal components of target genes") -> p

    print(p)
    dev.off()

    # visualise the betas somehow

    # Convert the list of matrices to a data frame for ggplot
    df <- do.call(rbind, lapply(seq_along(betas), function(x) {
      data.frame(matrix = x, as.data.frame(betas[[x]]))
    }))

    # Convert the data frame from wide to long format
    df %>% tidyr::pivot_longer(cols = starts_with("t")) -> df_long

    # Create a faceted histogram
    png(file.path(output_path, paste0("lm_data_gen_betas_hist.png")),
        width = 1024, height = 480, units = "px")

    ggplot(df_long, aes(x = value)) +
      geom_histogram(bins = 10
      ) +
      theme_minimal() +
      theme(#panel.grid.major = element_blank()#,
        panel.grid.minor = element_blank()
      )  +
      facet_grid(matrix ~ name, scales = "free", space = "free") +
      scale_x_continuous(breaks = c(-coefficients_standard_deviation,
                                    #0,
                                    coefficients_standard_deviation) ) +
      scale_y_continuous(breaks = NULL) +
      xlab("target genes") + ylab("cell clusters") -> p

    print(p)

    dev.off()


    # Visualise the mapping
    # since we generat everything as independent regressions we might as well
    # just plot one of the regression i guess
    vis_dat$residuals %>%
      tidyr::pivot_longer(cols = starts_with("t")) -> residuals_long

    # residuals_long  %>%
    #   group_by(name) %>%
    #   summarise(sd=sd(value), mean = mean(value)) -> temp_stat

    png(file.path(output_path, paste0("lm_data_gen_residuals.png")),
        width = 1024, height = 480, units = "px")

    ggplot(residuals_long, aes(x = value)) +
      geom_density(
      ) +
      theme_minimal() +
      theme(#panel.grid.major = element_blank()#,
        panel.grid.minor = element_blank()
      )  +
      facet_grid(cols = vars(name), scales = "free", space = "free") +
      scale_x_continuous(breaks = c(-target_gene_type_standard_deviation,
                                    #0,
                                    target_gene_type_standard_deviation) ) -> p
      print(p)
    dev.off()
  }

  return(list(dat = dat, true_betas = betas))

  # return indexes of different gene ?
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

  # Set seed for example
  set.seed(1234)

  dat <- generate_data_lm(n_cell_clusters = 3,
                          n_target_gene_type = 5,  # We have x named target genes that have one expression per cell
                          n_regulator_gene_type = 10,  # We have x named regulator genes that have one expression per cell
                          n_cells = c(1000, 2000, 30000),
                          regulator_means = pracma::linspace(1, 10, n = 3),  # Regulator mean expression in each cell cluster.
                          regulator_standard_deviations = pracma::linspace(0.1, 2, n = 3),  # Regulator sd for expression in each cell cluster.
                          coefficients_standard_deviation = 100, # 'The betas/slopes'. One per target gene. Instead of providing mean and standard deviation for each target gene, we provide the standard deviation from which these will be generated. Mean will be 0.
                          target_gene_type_standard_deviation = 3,
                          plot_stuff = T
  )
  true_betas <- dat$true_betas
  dat <- dat$dat

  dat[, 'true_cell_cluster_allocation'] <- paste("Cluster", pull(dat, 'true_cell_cluster_allocation'))  # These needs to be strings for discrete labels in pca plot
  pca_res <- prcomp(dat[, 3:ncol(dat)], scale. = TRUE)
  p <- ggplot2::autoplot(pca_res, data = dat, colour = 'true_cell_cluster_allocation')
  plot(p)
  print(dat)
  print(true_betas)
}



