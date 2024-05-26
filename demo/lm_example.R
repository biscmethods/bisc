#!/usr/bin/Rscript
rm(list = ls())

library(here)  # To work with paths
library(patchwork)

# options(warn = 2)  # To convert warning messages into error messages which display row of error. For debugging.

# Get absolute path where script is located, by using relative paths.
demo_path <- here::here("demo")
R_path <- here::here("R")
output_path <- demo_path

source(file.path(R_path, "generate_data_lm.R"))
source(file.path(R_path, "biclust_lm.R"))
source(file.path(R_path, "randomise_cluster_labels.R"))


# Set seed for example
set.seed(1234)

res <- generate_data_lm(n_cell_clusters = 3,
                        n_target_gene_type = 20,  # We have x named target genes that have one expression per cell
                        n_regulator_gene_type = 100,  # We have x named regulator genes that have one expression per cell
                        n_cells = c(1000, 2000, 3000),
                        regulator_means = c(1, 2, 3),  # Regulator mean expression in each cell cluster.
                        regulator_standard_deviations = c(0.1, 0.2, 0.2),  # Regulator sd for expression in each cell cluster.
                        coefficients_standard_deviation = 100, # 'The betas/slopes'. One per target gene. Instead of providing mean and standard deviation for each target gene, we provide the standard deviation from which these will be generated. Mean will be 0.
                        target_gene_type_standard_deviation = 3,
                        plot_stuff = TRUE
)

dat <- res$dat

# Split data into train/test
cell_data_split <- sample(c(1, 2), nrow(dat), replace = T)
train_indices <- which(cell_data_split == 1)
test_indices <- which(cell_data_split == 2)

dat_train <- dat[train_indices,]
dat_test <- dat[test_indices,]

disturbed_initial_cell_clust <- randomise_cluster_labels(cluster_labels = dat$true_cell_cluster_allocation,
                                                         fraction_randomised = 0.2)

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

penalization_lambdas <- sort(c(0, 10^(linspace(-6, 0, 10)), c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 2.0)))
rand_indexes_all <- vector(length = length(penalization_lambdas))
n_iterations_all <- vector(length = length(penalization_lambdas))
cluster_complexity_all <- vector(length = length(penalization_lambdas))
max_iter <- 15
BIC_all <- matrix(ncol = max_iter, nrow = length(penalization_lambdas))
target_genes_residual_var_all <- vector(mode = "list", length = length(penalization_lambdas))
for (i in seq_along(penalization_lambdas)) {
  penalization_lambda <- penalization_lambdas[i]
  penalization_lambda_str <- sprintf("%.6f", penalization_lambda)
  modded_output_path <- file.path(output_path, paste("penalization_lambda_", penalization_lambda_str))
  dir.create(modded_output_path, showWarnings = FALSE)

  res <- biclust(dat = dat_train,
                 dat_test = dat_test,
                 max_iter = max_iter,
                 initial_clustering = disturbed_initial_cell_clust_train,
                 n_target_genes = n_target_genes_train,
                 n_regulator_genes = n_regulator_genes_train,
                 n_total_cells = n_total_cells_train,
                 n_cell_clusters = n_cell_clusters,
                 ind_targetgenes = ind_targetgenes_train,
                 ind_reggenes = ind_reggenes_train,
                 output_path = modded_output_path,
                 penalization_lambda = penalization_lambda,
                 use_weights = TRUE)
  if (length(res) == 1 && is.na(res)) {
    rand_indexes_all[i] <- NA
    n_iterations_all[i] <- NA
    cluster_complexity_all[i] <- NA
    target_genes_residual_var_all[[i]] <- NA
    temp_BIC <- NA
    BIC_all[i, seq_along(temp_BIC)] <- NA
    print(paste("For penalization lambda:", penalization_lambda, ", algoritm crashes"), quote = FALSE)
  }else {
    rand_indexes_all[i] <- res$rand_index
    n_iterations_all[i] <- res$n_iterations
    cluster_complexity_all[i] <- res$db
    target_genes_residual_var_all[[i]] <- res$taget_genes_residual_var
    temp_BIC <- unlist(res$BIC)
    BIC_all[i, seq_along(temp_BIC)] <- temp_BIC
    print(paste("For penalization lambda:", penalization_lambda, ", Final rand index when compared to true clusters:", res$rand_index), quote = FALSE)
  }
  print("", quote = FALSE)

}

# Basic scatterplots of Penalization Lambdas vs Rand index
last_iteration_BIC <- unlist(apply(BIC_all, MARGIN = 1, FUN = function(x) x[max(which(!is.na(x)))][[1]], simplify = TRUE))

df <- data.frame("cluster_complexity" = as.numeric(as.character(cluster_complexity_all)),
                 "penalization_lambdas" = as.numeric(as.character(penalization_lambdas)),
                 "number_of_iterations" = as.numeric(as.character(n_iterations_all)),
                 "rand_index_result_vs_true" = as.numeric(as.character(rand_indexes_all)),
                 "BIC" = as.numeric(as.character(last_iteration_BIC))
)


df <- df[order(df$cluster_complexity, decreasing = TRUE),]

png(file.path(output_path, paste0("final_plot.png")),
    width = 1024, height = 480, units = "px")

p1 <- ggplot(data = df, aes(x = penalization_lambdas, y = cluster_complexity, group = 1)) +
  geom_line(color = "red") +
  geom_point() +
  scale_x_log10(trans = scales::pseudo_log_trans(sigma = 0.01)) +
  labs(x = "Penalization Lambda", y = "Silhoutte")

p2 <- ggplot(data = df, aes(x = number_of_iterations, y = cluster_complexity, group = 1)) +
  geom_line(color = "red") +
  geom_point() +
  labs(x = "Number of iterations", y = "Silhoutte")
p3 <- ggplot(data = df, aes(x = rand_index_result_vs_true, y = cluster_complexity, group = 1)) +
  geom_line(color = "red") +
  geom_point() +
  labs(x = "Rand Index, result vs true", y = "Silhoutte")
p4 <- ggplot(data = df, aes(x = penalization_lambdas, y = rand_index_result_vs_true, group = 1)) +
  geom_line(color = "red") +
  geom_point() +
  labs(x = "Penalization Lambda", y = "Rand Index, result vs true")
p5 <- ggplot(data = df, aes(x = penalization_lambdas, y = rand_index_result_vs_true, group = 1)) +
  geom_line(color = "red") +
  geom_point() +
  scale_x_log10(trans = scales::pseudo_log_trans(sigma = 0.01)) +
  labs(x = "Penalization Lambda", y = "Rand Index, result vs true")
p6 <- ggplot(data = df, aes(x = penalization_lambdas, y = BIC, group = 1)) +
  geom_line(color = "red") +
  geom_point() +
  scale_x_log10(trans = scales::pseudo_log_trans(sigma = 0.01)) +
  labs(x = "Penalization Lambda", y = "Last iteration BIC")
p1 + p2 + p3 + p4 + p5 + p6
dev.off()


