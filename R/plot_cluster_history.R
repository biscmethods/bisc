#!/usr/bin/Rscript

if (!require(aricode)) install.packages('aricode')
if (!require(ggplot2)) install.packages('ggplot2')
if (!require(dplyr)) install.packages('dplyr')
if (!require(ggalluvial)) install.packages('ggalluvial')
if (!require(reshape2)) install.packages('reshape2')
if (!require(ggsankey)) {
  if (!require(remotes)) install.packages('remotes')
  library(remotes)
  remotes::install_github("davidsjoberg/ggsankey")
}
library(ggsankey)
library(ggalluvial)
library(ggplot2)
library(reshape2)
library(dplyr)
library(aricode)  # To calculate rand index

plot_cluster_history <- function(cell_cluster_history, correct_plot = TRUE) {

  d <- cell_cluster_history
  d <- d[, colSums(is.na(d)) == 0]  # Remove NA data

  new_colnames <- colnames(d)

  rand_ind <- vector(length = (ncol(cell_cluster_history) - 1))
  for (i in 2:ncol(cell_cluster_history)) {
    rand_ind[i - 1] <- round(RI(cell_cluster_history[, 2],
                                cell_cluster_history[, i]), 2)
    new_colnames[i] <- paste0(new_colnames[i], "\nRI:", rand_ind[i - 1])
  }

  colnames(d) <- new_colnames

  d <- reshape2::melt(d, id.vars = "Cell ID")  # TODO: https://stackoverflow.com/questions/25688897/reshape2-melt-warning-message
  colnames(d) <- c("cell", "iteration", "cluster")
  d['cluster'] <- as.factor(d[, 'cluster'])

  # Plotting it
  # Slow. But keeps track of individual cells
  if (correct_plot) {
    p <- ggplot2::ggplot(d, ggplot2::aes(x = iteration,
                                         stratum = cluster,
                                         alluvium = cell,
                                         fill = cluster,
                                         label = cluster)) +
      ggplot2::scale_fill_brewer(type = "qual", palette = "Set2") +
      ggalluvial::geom_flow(stat = "alluvium", lode.guidance = "frontback") +
      ggalluvial::geom_stratum(alpha = 0.5) +
      ggplot2::geom_text(stat = "stratum", size = 3) +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      ggplot2::ggtitle("Alluvial diagram tracking cluster allocation in each iteration")
    print(p)
    # ggsave('cell_cluster_history.png', plot=p, dpi = 300, height = 6, width = 12, unit = 'in')
  }else {
    # Doesn't keep track of individual cells
    p <- ggplot2::ggplot(d, ggplot2::aes(x = iteration, stratum = cluster, alluvium = cell, fill = cluster, label = cluster)) +
      ggplot2::scale_fill_brewer(type = "qual", palette = "Set2") +
      ggalluvial::geom_flow() +
      ggalluvial::geom_stratum() +
      ggplot2::ylab("Cells") +
      ggplot2::xlab("Iteration") +
      ggplot2::labs(fill = "Cluster") +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::ggtitle(paste0("Log of cluster allocation\nRand index of true vs final: ", round(rand_ind, 2)))
    print(p)
  }


}
