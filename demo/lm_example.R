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
                        n_target_gene_type = 20,  # We have x named target genes that have one expression per cell
                        n_regulator_gene_type = 100,  # We have x named regulator genes that have one expression per cell
                        n_cells = c(1000, 2000, 3000),
                        regulator_means = c(1, 2, 3),  # Regulator mean expression in each cell cluster.
                        regulator_standard_deviations = c(0.1, 0.2, 0.2),  # Regulator sd for expression in each cell cluster.
                        coefficients_standard_deviation = 100, # 'The betas/slopes'. One per target gene. Instead of providing mean and standard deviation for each target gene, we provide the standard deviation from which these will be generated. Mean will be 0.
                        target_gene_type_standard_deviation = 3
)

# Split data into train/test
cell_data_split <- sample(c(1, 2), nrow(dat), replace = T)
train_indices <- which(cell_data_split == 1)
test_indices <- which(cell_data_split == 2)

dat_train <- dat[train_indices,]
dat_test <- dat[test_indices,]

disturbed_initial_cell_clust <- randomise_cluster_labels(cluster_labels = dat$true_cell_cluster_allocation,
                                                         fraction_randomised = 0.2)

print(str(disturbed_initial_cell_clust))

disturbed_initial_cell_clust_train <- disturbed_initial_cell_clust[train_indices]

# TODO: Put this inside generate_data_lm or something
# Setup variables that will be used throughout the script
# We assume target genes start with t then a number. And likewise for regulator genes.
n_total_cells_train <- nrow(dat_train)
ind_targetgenes_train <- which(str_detect(colnames(dat_train), "t\\d"))
ind_reggenes_train <- which(str_detect(colnames(dat_train), "r\\d"))

# Set up some variables
n_cell_clusters <- length(unique(disturbed_initial_cell_clust))
n_target_genes_train <- length(ind_targetgenes_train)
n_regulator_genes_train <- length(ind_reggenes_train)

penalization_lambdas <- c(0, 0.5, 1.0)
RI <- vector(length = length(penalization_lambdas))
for (i in seq_along(penalization_lambdas)) {
  penalization_lambda <- penalization_lambdas[i]
  penalization_lambda_str <- sprintf("%.2f", penalization_lambda)
  modded_output_path <- file.path(output_path, paste("penalization_lambda_", penalization_lambda_str))
  dir.create(modded_output_path, showWarnings = FALSE)

  ri <- biclust(dat = dat_train,
                max_iter = 50,
                initial_clustering = disturbed_initial_cell_clust_train,
                n_target_genes = n_target_genes_train,
                n_regulator_genes = n_regulator_genes_train,
                n_total_cells = n_total_cells_train,
                n_cell_clusters = n_cell_clusters,
                ind_targetgenes = ind_targetgenes_train,
                ind_reggenes = ind_reggenes_train,
                output_path = modded_output_path,
                penalization_lambda = penalization_lambda)
  RI[i] <- ri
  print(paste(penalization_lambda, ri))
}