#!/usr/bin/Rscript
rm(list = ls())

library(here)  # To work with paths
library(patchwork)
sink()

# options(warn=2)  # To convert warning messages into error messages which display row of error. For debugging.

# Get absolute path where script is located, by using relative paths.
demo_path <- here::here("demo")
R_path <- here::here("R")
path_data <- here::here('data')
output_path <- demo_path

path_env_data <- file.path(path_data, "env_data_pbmc10k_sctransform.RData")


source(file.path(R_path, "generate_dummy_data_for_cell_clustering.R"))
source(file.path(R_path, "bisc.R"))
source(file.path(R_path, "randomise_cluster_labels.R"))

# Set seed for example
set.seed(250)
load(path_env_data)

keep_cluster <- c(1, 2)
n_target_gene_modules <- c(5, 4)  # Number of target gene clusters in each cell cluster
n_cell_clusters <- length(keep_cluster)

logical_keep_cell_cluster <- true_cluster_allocation == 100
for(i_cell_cluster in keep_cluster){
  logical_keep_cell_cluster <- logical_keep_cell_cluster | (true_cluster_allocation==i_cell_cluster)
}

biclust_input_data <- d[,which(logical_keep_cell_cluster)]
true_cluster_allocation <- true_cluster_allocation[which(logical_keep_cell_cluster)]
true_cluster_allocation <- match(true_cluster_allocation, sort(unique(true_cluster_allocation)))


penalization_lambdas <- c(0.2)
# initial_clustering <- factor(sample(true_cluster_allocation))
initial_clustering <- factor(sample(unique(true_cluster_allocation), size=length(true_cluster_allocation), replace=TRUE))

sink("output_pbmc.txt", split=TRUE)


all_res <- vector(mode = "list", length = 1)
for (c_seed in seq(1)){
  set.seed(c_seed)
  BICLUST_RESULTS <- vector(mode = "list", length = length(penalization_lambdas))


  for (i_penalization_lambda in seq_along(penalization_lambdas)) {
    print("", quote = FALSE)
    print(paste("Running biclust for penalization_lambda", penalization_lambdas[i_penalization_lambda]), quote = FALSE)
    BICLUST_RESULTS[[i_penalization_lambda]] <- tryCatch(
      expr = bisc(dat = t(biclust_input_data),
                  cell_id = colnames(biclust_input_data),
                  true_cell_cluster_allocation = true_cluster_allocation,
                  max_iter = 10,
                  n_target_gene_clusters = n_target_gene_modules,
                  initial_clustering = initial_clustering,
                  n_cell_clusters = n_cell_clusters,
                  ind_targetgenes = ind_targetgenes,
                  ind_reggenes = ind_reggenes,
                  output_path = output_path,
                  penalization_lambda = penalization_lambdas[i_penalization_lambda],
                  use_complex_cluster_allocation = FALSE,
                  calculate_BIC = FALSE,
                  calculate_silhoutte = FALSE,
                  calculate_davies_bouldin_index = FALSE,
                  use_garbage_cluster_targets = FALSE),
      error = function(e) {
        warning(paste0("Error in bisc() for c_seed=", c_seed," lambda=", i_penalization_lambda, ": ", e$message))
        return(NULL)
      }
    )



  }
  all_res[[c_seed]] <- BICLUST_RESULTS

}
sink()

print("", quote = FALSE)
print("", quote = FALSE)
RIs <- rep(NA, length(all_res))
iterations <- rep(NA, length(all_res))
converged <- rep(FALSE, length(all_res))
isna <- rep(FALSE, length(all_res))
for(i_res in seq(length(all_res))){
  BICLUST_RESULTS <- all_res[[i_res]]
  for (i_penalization_lambda in seq_along(penalization_lambdas)) {
    if (is.na(BICLUST_RESULTS[i_penalization_lambda])) {
      cat(paste("NA\n"))
      RIs[i_res] <- NA
      isna[i_res] <- TRUE
    } else if (is.null(BICLUST_RESULTS[i_penalization_lambda])) {
      cat(paste("NULL\n"))
      RIs[i_res] <- NULL
      isna[i_res] <- TRUE
    }else {

      RIs[i_res] <- BICLUST_RESULTS[[i_penalization_lambda]]$rand_index
      iterations[i_res] <- BICLUST_RESULTS[[i_penalization_lambda]]$n_iterations
      converged[i_res] <- BICLUST_RESULTS[[i_penalization_lambda]]$converged

      cat(paste0(BICLUST_RESULTS[[i_penalization_lambda]]$rand_index, ",", iterations[i_res], ",", converged[i_res]))
      cat("\n")

    }
  }
}


nonconverged_RIs <- RIs[!converged & !isna]
mean(nonconverged_RIs)
converged_RIs <- RIs[converged]
mean(converged_RIs)

svg(file.path(output_path, paste0("pbmc10k_RI.svg")), width=8, height=6)
boxplot(RIs ~ converged,
        ylab = "RIs",
        xlab = "Convergence Status",
        main = "Distribution of rand indexes (RIs)\nby whether or not they converged or not.",
        col = c("lightblue", "lightpink"))
dev.off()

svg(file.path(output_path, paste0("pbmc10k_iterations.svg")), width=8, height=6)
boxplot(iterations,
        ylab = "Iterations",
        xlab = "",
        main = "Iterations in bisc needed to converge.",
        col = c("lightblue", "lightpink"))
dev.off()
