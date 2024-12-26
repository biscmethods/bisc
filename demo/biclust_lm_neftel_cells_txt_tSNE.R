#!/usr/bin/Rscript

rm(list = ls())

library(here)  # To work with paths
library(patchwork)
library(tidyverse)
# sink()

options(warn=2)  # To convert warning messages into error messages which display row of error. For debugging.

# Get absolute path where script is located, by using relative paths.
demo_path <- here::here("demo")
R_path <- here::here("R")
output_path <- demo_path

source(file.path(R_path, "biclust_lm.R"))
source(file.path(R_path, "randomise_cluster_labels.R"))

#############################################
############ load data    ###################
#############################################

# Set seed for example
set.seed(214)

# Load data
# Folders and filenames
path_data <- here::here('data')
path_Neftel2019 <- file.path(path_data, "Neftel2019")
path_Group2 <- file.path(path_Neftel2019, "Group2")

cells <- read.table(file = file.path(path_Group2, 'cells.csv'),
                    sep = ' ',
                    header = TRUE,
                    stringsAsFactors = FALSE) %>% as_tibble()

# sanity check
for(iter in seq_along(unique(cells$cell_type))[c(1,3)]
    ){
  type <- unique(cells$cell_type)[iter]
  print(type)
  cells %>% filter(cell_type == type) %>% lm(tSNE1 ~ tSNE2,
                                           data = .) -> mm
  print(mm)
}

# use these cell types, others are meh
working_cell_types <- unique(cells$cell_type)[c(1,3)]

true_cell_cluster_allocation <- vector(length= length(cells$MESlike2)) * NA

for(iter in seq_along(unique(
  cells$cell_type[(cells$cell_type %in% working_cell_types)]))
    )
{
  type <- unique(
    cells$cell_type[(cells$cell_type %in% working_cell_types)])[iter]
  true_cell_cluster_allocation[cells$cell_type == type] <- iter
}

unique(true_cell_cluster_allocation[(cells$cell_type %in% working_cell_types)])


dat <- tibble (cell_id = cells$cell_name,
               true_cell_cluster_allocation,
               cells$tSNE1,
               cells$tSNE2
               #, malignant = (cells$malignant == "yes") + 0
               ) %>% filter(
                 cell_id %in%
                   cells$cell_name[
                     which(cells$cell_type %in% working_cell_types)
                     ]
                 )
# dat %>% filter(!is.na(malignant)) -> dat


# Split data into train/test
cell_data_split <- sample(c(1, 2), nrow(dat), replace = T)
train_indices <- which(cell_data_split == 1)
test_indices <- which(cell_data_split == 2)

dat_train <- dat[train_indices,]
dat_test <- dat[test_indices,]

disturbed_initial_cell_clust <- randomise_cluster_labels(cluster_labels = dat$true_cell_cluster_allocation,
                                                         fraction_randomised = 0.5)

disturbed_initial_cell_clust_train <- disturbed_initial_cell_clust[train_indices]

# Setup variables that will be used throughout the script
# We assume target genes start with t then a number. And likewise for regulator genes.
n_total_cells_train   <- nrow(dat_train)
ind_targetgenes_train <- which(str_detect(colnames(dat_train), "tSNE2"))
ind_reggenes_train    <- which(str_detect(colnames(dat_train), "tSNE1"))

# n_total_cells_test    <- nrow(dat_test)

# Set up some variables
n_cell_clusters <- length(unique(disturbed_initial_cell_clust))
n_target_genes_train <- length(ind_targetgenes_train)
n_regulator_genes_train <- length(ind_reggenes_train)

penalization_lambdas <- sort(c( 10^(linspace(-6, 0, 10)), c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 2.0)))
rand_indexes_all <- vector(length = length(penalization_lambdas))
n_iterations_all <- vector(length = length(penalization_lambdas))
cluster_complexity_all <- vector(length = length(penalization_lambdas))
max_iter <- 15
BIC_all <- matrix(ncol = max_iter, nrow = length(penalization_lambdas))
target_genes_residual_var_all <- vector(mode = "list", length = length(penalization_lambdas))


for (i in seq_along(penalization_lambdas)) {
  penalization_lambda <- penalization_lambdas[i]
  penalization_lambda_str <- sprintf("%.6f", penalization_lambda)
  modded_output_path <- file.path(output_path, paste("biclust_lm_neftel_penalization_lambda_", penalization_lambda_str))
  dir.create(modded_output_path, showWarnings = FALSE)

  tryCatch({
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
                 use_weights = FALSE,
                 use_complex_cluster_allocation = FALSE
                 )
  },
  error = function(e) {
    print("An error has occurred:")
    print(e)
    print("Traceback:")
    traceback()
  })



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

df <- data.frame(
                 "cluster_complexity" = as.numeric(as.character(cluster_complexity_all)),
                 "penalization_lambdas" = as.numeric(as.character(penalization_lambdas)),
                 "number_of_iterations" = as.numeric(as.character(n_iterations_all)),
                 "rand_index_result_vs_true" = as.numeric(as.character(rand_indexes_all)),
                 "BIC" = as.numeric(as.character(last_iteration_BIC))
)


df <- df[order(df$cluster_complexity, decreasing = TRUE),]

png(file.path(output_path, paste0("final_plot_lm_neftel_cells_txt.png")),
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



# Runs only when script is run by itself
# || interactive()
if (sys.nframe() == 0) {
  p1 + p2 + p3 + p4 + p5 + p6
}
