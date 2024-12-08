#!/usr/bin/Rscript

library(here)  # To work with paths
library(patchwork)

while (sink.number() > 0) {
  sink()
  sink(file = NULL)
}
gc()

# options(warn=2)  # To convert warning messages into error messages which display row of error. For debugging.

raw_printoutput_path <- file.path(local_data, "output_pbmc_300.txt")
all_res_path <- file.path(local_data, "run300.rds")

source(file.path(R_path, "generate_dummy_data_for_cell_clustering.R"))
source(file.path(R_path, "bisc.R"))
source(file.path(R_path, "randomise_cluster_labels.R"))

# Set seed for example
set.seed(250)

n_target_gene_modules <- c(5, 4)  # Number of target gene clusters in each cell cluster

if(exists('true_cluster_allocation')){
  n_cell_clusters <- length(unique(true_cluster_allocation))
}else{
  n_cell_clusters <- length(n_target_gene_modules)
}

seeds <- seq(300)

penalization_lambdas <- c(0.2)

if(!file.exists(raw_printoutput_path) || !file.exists(all_res_path) || redo_flag){

  # initial_clustering <- factor(sample(true_cluster_allocation))
  initial_clustering <- factor(sample(unique(true_cluster_allocation), size=length(true_cluster_allocation), replace=TRUE))



  sink(raw_printoutput_path, split=TRUE)

  all_res <- vector(mode = "list", length = length(seeds))
  for (i_seed in seq_along(seeds)){
    current_seed <- seeds[i_seed]
    set.seed(current_seed)
    BICLUST_RESULTS <- vector(mode = "list", length = length(penalization_lambdas))

    for (i_penalization_lambda in seq_along(penalization_lambdas)) {
      print("", quote = FALSE)
      print(paste("Running biclust for penalization_lambda", penalization_lambdas[i_penalization_lambda]), quote = FALSE)
      BICLUST_RESULTS[[i_penalization_lambda]] <- tryCatch(
        expr = bisc(dat = t(d),
                    cell_id = colnames(d),
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
                    calculate_BIC = TRUE,
                    calculate_silhoutte = TRUE,
                    calculate_davies_bouldin_index = TRUE,
                    use_garbage_cluster_targets = FALSE),
        error = function(e) {
          warning(paste0("Error in bisc() for c_seed=", current_seed," lambda=", i_penalization_lambda, ": ", e$message))
          return(NULL)
        }
      )



    }
    all_res[[i_seed]] <- BICLUST_RESULTS

  }
  saveRDS(all_res, all_res_path)
  sink()
}else{
  all_res <- readRDS(all_res_path)
}


# Extract data from all_res ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("\n\nExtract data\n")
RIs <- rep(NA, length(all_res))
iterations <- rep(NA, length(all_res))
converged <- rep(FALSE, length(all_res))
isna <- rep(FALSE, length(all_res))
for(i_res in seq(length(all_res))){
  current_seed <- seeds[i_res]
  BICLUST_RESULTS <- all_res[[i_res]]
  for (i_penalization_lambda in seq_along(penalization_lambdas)) {
    if (is.na(BICLUST_RESULTS[i_penalization_lambda])) {
      cat(paste0(" Seed:", current_seed, " is NA\n"))
      RIs[i_res] <- NA
      isna[i_res] <- TRUE
    } else if (is.null(BICLUST_RESULTS[i_penalization_lambda])) {
      cat(paste0(" Seed:", current_seed, " is NULL\n"))
      RIs[i_res] <- NULL
      isna[i_res] <- TRUE
    }else {

      RIs[i_res] <- BICLUST_RESULTS[[i_penalization_lambda]]$rand_index
      iterations[i_res] <- BICLUST_RESULTS[[i_penalization_lambda]]$n_iterations
      converged[i_res] <- BICLUST_RESULTS[[i_penalization_lambda]]$converged

      cat(paste0(" Seed:", current_seed, " is RI:", BICLUST_RESULTS[[i_penalization_lambda]]$rand_index, ", n_iter:", iterations[i_res], ", converged:", converged[i_res]))
      cat("\n")

    }
  }
}


nonconverged_RIs <- RIs[!converged & !isna]
# mean(nonconverged_RIs)
converged_RIs <- RIs[converged]
# mean(converged_RIs)


# Plot info -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
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



# plot_heatmap ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cell_clustering_RIs <- RIs[converged]
min(RIs[converged])
max(RIs[converged])
current_biclust <- (all_res[[which(RIs==max(RIs, na.rm=TRUE))]][[1]])
target_gene_allocation <- current_biclust$scregclust_final_result_module
cell_cluster_allocation <- current_biclust$cell_cluster_allocation
# Assign one unique number to each gene module --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if(exists('n_total_cells') & exists('n_target_genes') ){
  biclust_result_matrix <- matrix(0, nrow = n_total_cells, ncol = n_target_genes)
}else{
  n_total_cells <- length(all_res[[2]][[1]]$call$initial_clustering)

  n_target_genes <- length(all_res[[2]][[1]]$call$ind_targetgenes)

  biclust_result_matrix <- matrix(0, nrow = n_total_cells, ncol = n_target_genes)
}


# Check each entry in the matrix (each cell and gene pair),
# and assign a string-number to each unique "cell-cluster-gene-module".
for (i in 1:n_total_cells) {

  if(is.na(cell_cluster_allocation[i])){
    cat("Cell cluster allocation is NA. Random cluster will be assigned.\n")
    cell_cluster_allocation[i] <- sample.int(length(target_gene_allocation),1)
  }

  cluster <- cell_cluster_allocation[i]
  gene_allocation <- target_gene_allocation[[cluster]][1:n_target_genes]
  gene_allocation[gene_allocation==-1] <- 0
  # Create unique numbers for each pair
  biclust_result_matrix[i, ] <- paste0(cluster, gene_allocation)
}

# Convert the string-numbers to numeric matrix (starting from 1 this time)
biclust_result_matrix <- matrix(as.numeric(biclust_result_matrix),
                                nrow = n_total_cells,
                                ncol = n_target_genes)

biclust_result_matrix <- matrix(as.integer(as.factor(biclust_result_matrix)),
                                nrow=nrow(biclust_result_matrix),
                                ncol=ncol(biclust_result_matrix))


# Make the plots --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

n <- length(unique(as.vector(biclust_result_matrix)))
regions <- seq(1, n, length.out = n + 1)
middle_of_regions <- (regions[-1] + regions[-length(regions)]) / 2
odd_number_larger <- ifelse(n %% 2 == 0, n + 1, n)
if(n %% 2 == 1){
  keep_these_colors = 1:n
}else{
  keep_these_colors <- setdiff(1:odd_number_larger, (odd_number_larger + 1) / 2 + 1)
}
constructed_plot <- rasterVis::levelplot(biclust_result_matrix,
                                       att = n,
                                       # col.regions = rainbow(odd_number_larger),
                                       colorkey = list(at = regions,
                                                       # col=rainbow(odd_number_larger)[keep_these_colors],
                                                       labels = list(at = middle_of_regions, labels = as.character(1:n))),
                                       xlab = 'Cells',
                                       ylab = 'Target genes')


# Plot in IDE
print(constructed_plot)
aspect_ratio <- nrow(biclust_result_matrix) / ncol(biclust_result_matrix)
print(aspect_ratio)
# Plot to file
png(file.path(output_path, paste0("pbmc_bisc_heatmap.png")), width = 100*aspect_ratio, height = 200, units = "px")
print(constructed_plot)
dev.off()
pdf(file.path(output_path, paste0("pbmc_bisc_heatmap.pdf")), width =  0.95*aspect_ratio/2, height = 1.9)
print(constructed_plot)
dev.off()
