#!/usr/bin/Rscript

library(tidyverse)
library(aricode)  # To calculate rand index
library(ggplot2)
library(ggalluvial)
library(reshape2)
library(ggfortify)  # For pca-plot
library(stats)  # For prcomp


# Scatter-plot the log-liklihood on each axis, color with true allocation
# If more than 2 dim make pca plot
scatter_plot_loglikelihood <- function(dat,  # cells x (target_genes + regulator_genes)
                                       likelihood,  # cells x cell_clusters
                                       n_cell_clusters,
                                       penalization_lambda,
                                       output_path,
                                       i_main,
                                       true_cell_cluster_allocation_vector = NULL) {

  if (is.null(true_cell_cluster_allocation_vector) && ('true_cell_cluster_allocation' %in% colnames(dat))) {
    true_cell_cluster_allocation_vector <- paste("Cluster", pull(dat, var = 'true_cell_cluster_allocation'))  # These needs to be strings for discrete labels in pca plot
  }else if (is.null(true_cell_cluster_allocation_vector)) {
    stop("True cell clusters needs to be supplied as a column in dat or as a seperate variable")
  }

  if (is.matrix(likelihood)) {
    colnames(likelihood) <- paste("Likelihood cell cluster", seq_len(n_cell_clusters))
  }
  likelihood_tibble <- tibble::as_tibble(likelihood)
  data_for_plotting <- tibble::tibble(likelihood_tibble, true_cell_cluster_allocation = true_cell_cluster_allocation_vector)

  filename_plot <- paste0("Decision_line_lambda_", round(penalization_lambda, digits = 3), "_iteration_", i_main, ".png")
  if (!is.matrix(likelihood) || ncol(likelihood) == 0) {
    p <- ggplot2::ggplot(data = data_for_plotting, ggplot2::aes(x = "Likelihood cell cluster 1", y = "Likelihood cell cluster 2", color = true_cell_cluster_allocation)) +
      ggplot2::geom_point() +
      ggplot2::geom_abline(intercept = 0, slope = 1) +
      ggplot2::labs(x = "Log-likelihood for fitting into cluster 1", y = "Log-likelihood for fitting into cluster 2") +
      ggplot2::labs(color = "True cell cluster")
    png(file.path(output_path, filename_plot))
    plot(p)
    dev.off()
  } else {
    if (any(sapply(data_for_plotting[, seq_len(ncol(data_for_plotting) - 1)], function(x) all(x == 0)))) {
      pca_res <- stats::prcomp(data_for_plotting[, seq_len(ncol(data_for_plotting) - 1)], scale. = FALSE)
    }else {
      pca_res <- stats::prcomp(data_for_plotting[, seq_len(ncol(data_for_plotting) - 1)], scale. = TRUE)
    }
    p <- ggplot2::autoplot(pca_res, data = data_for_plotting, colour = 'true_cell_cluster_allocation', alpha = 0.1)
    png(file.path(output_path, filename_plot))
    plot(p)
    dev.off()
  }
}


# Make histograms
hist_plot_loglikelihood <- function(dat,
                                    likelihood,
                                    n_cell_clusters,
                                    penalization_lambda,
                                    output_path,
                                    i_main,
                                    true_cell_cluster_allocation_vector = NULL) {
  if (is.null(true_cell_cluster_allocation_vector) && ('true_cell_cluster_allocation' %in% colnames(dat))) {
    true_cell_cluster_allocation_vector <- pull(dat, var = 'true_cell_cluster_allocation')
  }else if (is.null(true_cell_cluster_allocation_vector)) {
    stop("True cell clusters needs to be supplied as a column in dat or as a seperate variable")
  }

  colnames(likelihood) <- paste("Likelihood cell cluster", seq_len(n_cell_clusters))
  likelihood_tibble <- tibble::as_tibble(likelihood)

  likelihood_tibble['cell_id'] <- seq_len(nrow(likelihood_tibble))
  plots <- vector(mode = "list", length = n_cell_clusters)
  for (cell_cluster in seq_len(n_cell_clusters)) {
    cell_cluster_rows <- which(true_cell_cluster_allocation_vector == cell_cluster)

    cell_cluster_likelihood <- likelihood_tibble[cell_cluster_rows,]
    # data_for_plotting <- tibble::tibble(cell_cluster_likelihood, true_cell_cluster_allocation = true_cell_cluster_allocation_vector)
    cell_cluster_likelihood <- reshape2::melt(cell_cluster_likelihood, id.vars = "cell_id")

    # Interleaved histograms
    plots[[cell_cluster]] <- ggplot2::ggplot(cell_cluster_likelihood, ggplot2::aes(x = value, color = variable)) +
      ggplot2::geom_histogram(fill = "white", position = "dodge", bins = 100) +
      ggplot2::theme(legend.position = "top") +
      ggplot2::ggtitle(label = paste("True cell cluster", cell_cluster)) #+
    #ggplot2::coord_trans(x = "log2")
  }
  title <- cowplot::ggdraw() + cowplot::draw_label("Likelihoods for the cells belonging in each true cell cluster.", fontface = 'bold')
  p <- cowplot::plot_grid(
    plotlist = plots,
    align = "hv"
  )
  p <- cowplot::plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control title margins
  filename_plot <- paste0("Histograms_lambda_", round(penalization_lambda, digits = 3), "_iteration_", i_main, ".png")
  png(file.path(output_path, filename_plot), width = 1024, height = 800)
  plot(p)
  dev.off()
}
