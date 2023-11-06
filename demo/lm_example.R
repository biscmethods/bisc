#!/usr/bin/Rscript
rm(list = ls())

library(here)  # To work with paths
# options(warn=2)  # To convert warning messages into error messages which display row of error. For debugging.

# Get absolute path where script is located, by using relative paths.
demo_path <- here::here("demo")
R_path <- here::here("R")
output_path <- demo_path

source(file.path(R_path, "generate_data_lm.R"))
source(file.path(R_path, "biclust.R"))
source(file.path(R_path, "randomise_cluster_labels.R"))

# Set seed for example
set.seed(1234)

dat <- generate_data_lm(n_cell_clusters = 3,
                        n_target_gene_type = 4,  # We have x named target genes that have one expression per cell
                        n_regulator_gene_type = 20,  # We have x named regulator genes that have one expression per cell
                        n_cells = c(1000, 5000, 4000),
                        regulator_means = c(1, 2, 3),  # Regulator mean expression in each cell cluster.
                        regulator_standard_deviations = c(0.1, 0.2, 0.2),  # Regulator sd for expression in each cell cluster.
                        coefficients_standard_deviation = 100, # 'The betas/slopes'. One per target gene. Instead of providing mean and standard deviation for each target gene, we provide the standard deviation from which these will be generated. Mean will be 0.
                        target_gene_type_standard_deviation = 3
)

# TODO: Put this inside generate_data_lm or something
# Setup variables that will be used throughout the script
# We assume target genes start with t then a number. And likewise for regulator genes.
n_total_cells <- nrow(dat)
ind_targetgenes <- which(str_detect(colnames(dat), "t\\d"))
ind_reggenes <- which(str_detect(colnames(dat), "r\\d"))

disturbed_initial_cell_clust <- randomise_cluster_labels(cluster_labels = dat$true_cell_cluster_allocation)

# Set up some variables
n_cell_clusters <- length(unique(disturbed_initial_cell_clust))
n_target_genes <- length(ind_targetgenes)
n_regulator_genes <- length(ind_reggenes)

biclust(dat = dat,
        max_iter = 50,
        initial_clustering = disturbed_initial_cell_clust,
        n_target_genes = n_target_genes,
        n_regulator_genes = n_regulator_genes,
        n_total_cells = n_total_cells,
        n_cell_clusters = n_cell_clusters,
        ind_targetgenes = ind_targetgenes,
        ind_reggenes = ind_reggenes,
        output_path = output_path,
        penalization_lambda = 0.5)