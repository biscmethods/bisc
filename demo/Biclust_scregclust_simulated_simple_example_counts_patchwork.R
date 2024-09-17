#!/usr/bin/Rscript
rm(list = ls())

library(here)  # To work with paths
library(patchwork)
library(biclust)
library(rasterVis)
library(cowplot)
sink()

# options(warn=2)  # To convert warning messages into error messages which display row of error. For debugging.

# Get absolute path where script is located, by using relative paths.
demo_path <- here::here("demo")
R_path <- here::here("R")
output_path <- demo_path
path_data <- here::here('data')

redo_flag <- TRUE

source(file.path(R_path, "generate_dummy_data_for_cell_clustering.R"))
source(file.path(R_path, "biclust_scregclust.R"))


#############################################
############ data for dev ###################
#############################################

# Set seed for example
set.seed(1234)

# Set variables ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
n_cell_clusters <- 2
n_target_gene_clusters <- c(4, 2)  # Number of target gene clusters in each cell cluster
n_target_genes <- 100
n_regulator_genes <- 10
n_cells <- c(1000, 1000)
regulator_means <- c(0, 0) # For generating dummy data, regulator mean in each cell cluster
coefficient_means <- list(c(1, 3, 5, 7), c(10, 20))  # For generating dummy data, coefficient means in each cell cluster
coefficient_sds <- list(c(0.01, 0.01, 0.01, 0.01), c(0.01, 0.01))
disturbed_fraction <- 0.1  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
testing_penalization_data_gen <- c(0.1, 0.5)

if (!file.exists(file.path(path_data, "env_sim_simple_data_biclust_sc.rds")) |
    redo_flag) {
  generated_data <- generate_dummy_data_for_cell_clustering(
    n_cell_clusters = n_cell_clusters,
    n_target_gene_clusters = n_target_gene_clusters,
    # Number of target gene clusters in each cell cluster
    n_target_genes = n_target_genes,
    #from vignette
    n_regulator_genes = n_regulator_genes,
    # from vignette
    n_cells = n_cells,
    regulator_means = regulator_means,
    # For generating dummy data, regulator mean in each cell cluster
    coefficient_means <- coefficient_means,
    # For generating dummy data, coefficient means in each cell cluster
    coefficient_sds <- coefficient_sds,
    disturbed_fraction = disturbed_fraction,
    # Value between 0 and 1. How large portion of cells should move to other cell clusters.
    plot_stuff = FALSE,
    plot_suffix = "Simple",
    testing_penalization = testing_penalization_data_gen
  )


  saveRDS(generated_data,
          file.path(path_data, "env_sim_simple_data_biclust_sc.rds"))

} else {
  generated_data <- readRDS(file.path(path_data, "env_sim_simple_data_biclust_sc.rds"))

}

# plot counts
# devtools::install_github("jokergoo/ComplexHeatmap")
d <- generated_data$counts
d_u <- sort(unique(as.vector(d)))
colors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlBu"))(max(length(d_u), max(d_u)))
colors[1] <- "#000000"
colors <- structure(colors, names = as.character(d_u))
ComplexHeatmap::Heatmap(d,
                        col = colors,
                        name = "Generated data counts",
                        cluster_rows = FALSE,
                        cluster_columns = FALSE,
                        column_title = paste(ncol(d), "Genes"),
                        row_title = paste(nrow(d), "Cells"))


# library(Seurat)
# library(sctransform)
# Cells as columns according to:
# https://satijalab.github.io/seurat-object/reference/CreateSeuratObject.html
simdata <- Seurat::CreateSeuratObject(counts = t(generated_data$counts))
simdata <- Seurat::SCTransform(simdata)
simdata@assays$SCT@SCTModel.list$counts@feature.attributes$theta
simdata <- simdata@assays$SCT$scale.data



#------------------------------------------------------------------

# Because "dat <- cbind(Z_t, Z_r)" in generate_dummy_data_for_cell_clustering
ind_targetgenes <- which(c(rep(1, n_target_genes), rep(0, n_regulator_genes)) == 1)
ind_reggenes <- which(c(rep(0, n_target_genes), rep(1, n_regulator_genes)) == 1)


disturbed_initial_cell_clust <- factor(generated_data$disturbed_initial_cell_clust)

biclust_input_data <- t(simdata)
colnames(biclust_input_data) <- c(paste0("t", 1:n_target_genes),
                                  paste0("r", 1:n_regulator_genes))
biclust_input_data <- tibble::as_tibble(biclust_input_data)

# Set up some variables
n_cell_clusters <- length(unique(disturbed_initial_cell_clust))
n_target_genes <- length(ind_targetgenes)
n_regulator_genes <- length(ind_reggenes)
n_total_cells <- sum(n_cells)
cell_id <- 1:n_total_cells

#############################################
############ end data for dev ###############
#############################################


penalization_lambdas <- c(0.0001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7) # c( 0.00001, 0.1, 0.2, 0.5)
BICLUST_RESULTS <- vector(mode = "list", length = length(penalization_lambdas))

if (!file.exists(file.path(path_data, "env_sim_simple_nogarb_res_biclust_sc.rds")) |
    redo_flag) {
  for (i_penalization_lambda in seq_along(penalization_lambdas)) {
    print("", quote = FALSE)
    print(paste(
      "Running biclust for penalization_lambda",
      penalization_lambdas[i_penalization_lambda]
    ),
    quote = FALSE)
    BICLUST_RESULTS[[i_penalization_lambda]] <- biclust_scregclust(
      dat = biclust_input_data,
      cell_id = cell_id,
      true_cell_cluster_allocation = factor(generated_data$true_cell_clust),
      max_iter = 100,
      n_target_gene_clusters = n_target_gene_clusters,
      initial_clustering = disturbed_initial_cell_clust,
      n_cell_clusters = n_cell_clusters,
      ind_targetgenes = ind_targetgenes,
      ind_reggenes = ind_reggenes,
      output_path = output_path,
      penalization_lambda = penalization_lambdas[i_penalization_lambda],
      use_complex_cluster_allocation = FALSE,
      calculate_BIC = FALSE,
      calculate_silhoutte = FALSE,
      calculate_davies_bouldin_index = FALSE,
      plot_suffix = "Simple_cluster_all",
      always_use_flat_prior = FALSE,
      use_garbage_cluster_targets = FALSE
    )
  }

  saveRDS(
    BICLUST_RESULTS,
    file.path(path_data, "env_sim_simple_nogarb_res_biclust_sc.rds")
  )

} else {
  BICLUST_RESULTS <- readRDS(file.path(path_data, "env_sim_simple_nogarb_res_biclust_sc.rds"))

}

print("", quote = FALSE)
print("", quote = FALSE)
for (i_penalization_lambda in seq_along(penalization_lambdas)) {
  if (is.na(BICLUST_RESULTS[i_penalization_lambda])) {
    print(paste("penalization_lambda", penalization_lambdas[i_penalization_lambda], "is NA"),
          quote = FALSE)
  } else if (is.null(BICLUST_RESULTS[i_penalization_lambda])) {
    print(paste("penalization_lambda", penalization_lambdas[i_penalization_lambda], "is NULL"),
          quote = FALSE)
  } else {
    print(
      paste(
        "penalization_lambda",
        penalization_lambdas[i_penalization_lambda],
        "is ok with rand index",
        BICLUST_RESULTS[[i_penalization_lambda]]$rand_index
      ),
      quote = FALSE
    )
  }
}


##############################
# run again many times for stats
###########################

# compare with other biclust variant
standard_biclust_results <- biclust::biclust(
  x = as.matrix(biclust_input_data),
  method = BCSpectral(),
  normalization = "irrc",
  numberOfEigenvalues = 6,
  minr = 2,
  minc = 2,
  withinVar = 1
)


res1 <- biclust::biclust(
  as.matrix(biclust_input_data[, 1:n_target_genes]),
  method = BCCC(),
  delta = 1.5,
  alpha = 1,
  number = 2
)
res1 <- biclust::biclust(
  as.matrix(biclust_input_data[, 1:n_target_genes]),
  method=BCPlaid(),
  background=FALSE,
  iter.startup=100,
  iter.layer=100,
  back.fit=100,
  row.release=0.7,
  col.release=0.7,
  shuffle=100,
  max.layers=5,
  verbose=FALSE)

# a <- as.matrix(biclust_input_data)
# # colnames(a) <- 1:ncol(a)
# png("heatmap.png",
#     width = 8000,
#     height = 7000,
#     res = 300)
# biclust::heatmapBC(
#   x = a,
#   bicResult = res1,
#   order = FALSE,
#   Rowv = FALSE,
#   Colv = FALSE,
#   labRow = NA,
#   labCol = NA
# )
# dev.off()

calc_hamming <- function(matrix_data){
  hamming_dist <- function(vec1, vec2) {
    sum(vec1 != vec2)
  }

  # Number of rows (vectors)
  n <- nrow(matrix_data)

  # Initialize a distance matrix
  dist_matrix <- matrix(0, n, n)

  # Calculate pairwise Hamming distances between rows
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      dist_matrix[i, j] <- hamming_dist(matrix_data[i, ], matrix_data[j, ])
      dist_matrix[j, i] <- dist_matrix[i, j]  # Symmetric matrix
    }
  }

  # Print the distance matrix
  print(dist_matrix)
}

b <- as.matrix(biclust_input_data)[,1:n_target_genes]
b <- matrix(0, nrow = nrow(b), ncol = ncol(b))
for (i_n in 1:res1@Number) {
  a <- as.matrix((res1@RowxNumber[, i_n, drop = FALSE] %*% res1@NumberxCol[i_n, , drop = FALSE]) == 1)
  b[a] <- i_n
}
matrix_data <- unique(b, MARGIN = 1)
if(nrow(matrix_data)==1){
  print("only one type of vector")
}else{
  dist_matrix <- calc_hamming(matrix_data)
  distance_object <- as.dist(dist_matrix)

  # Perform hierarchical clustering
  hc <- hclust(distance_object)

  # Plot the dendrogram to visualize the clustering
  plot(hc, labels = rownames(matrix_data))
  clusters <- cutree(hc, k = 3)
  clusters[clusters>2] = 2

  # Print the cluster assignments
  print(clusters)

  res_cell_cluster <- vector(length=nrow(b))
  for (i in 1:nrow(matrix_data)) {
    inds <- which(apply(b, 1, function(x) return(all(x == matrix_data[i,]))))
    res_cell_cluster[inds] <- clusters[i]
  }



  # b <- raster::ratify(raster::raster(b))
  n <- length(unique(as.vector(b)))
  regions <- seq(1, n, length.out = n + 1)
  middle_of_regions <- (regions[-1] + regions[-length(regions)]) / 2
  odd_number_larger <- ifelse(n %% 2 == 0, n + 1, n)
  if(n %% 2 == 1){
    keep_these_colors = 1:n
  }else{
    keep_these_colors <- setdiff(1:odd_number_larger, (odd_number_larger + 1) / 2 + 1)
  }
  rasterVis::levelplot(b, att = n,
                       # col.regions = rainbow(odd_number_larger),
                       colorkey = list(at = regions,
                                       # col=rainbow(odd_number_larger)[keep_these_colors],
                                       labels = list(at = middle_of_regions, labels = as.character(1:n))),
                       xlab = 'Cells',
                       ylab = 'Target genes',
                       main='biclust')
}
#----- Construct heatmap for generated data
res <- generated_data
cell_cluster_allocation <- res$true_cell_clust
target_gene_allocation <- res$true_target_gene_allocation

result_matrix <- matrix(0, nrow = n_total_cells, ncol = n_target_genes)

for (i in 1:n_total_cells) {
  cluster <- cell_cluster_allocation[i]
  gene_allocation <- target_gene_allocation[[cluster]][1:n_target_genes]
  gene_allocation[gene_allocation==-1] <- 0
  # Create unique numbers for each pair
  result_matrix[i, ] <- paste0(cluster, gene_allocation)
}

# Convert the result to numeric matrix
result_matrix <- matrix(as.numeric(result_matrix), nrow = n_total_cells, ncol = n_target_genes)
result_matrix <- matrix(as.integer(as.factor(result_matrix)), nrow=nrow(result_matrix), ncol=ncol(result_matrix))
correct_clustering <- as.vector(result_matrix)
calc_hamming(unique(result_matrix, MARGIN = 1))


n <- length(unique(as.vector(result_matrix)))
if(n!=max(as.vector(result_matrix))){
  print("Warning")
}
regions <- seq(1, n, length.out = n + 1)
middle_of_regions <- (regions[-1] + regions[-length(regions)]) / 2
odd_number_larger <- ifelse(n %% 2 == 0, n + 1, n)
if(n %% 2 == 1){
  keep_these_colors = 1:n
}else{
  keep_these_colors <- setdiff(1:odd_number_larger, (odd_number_larger + 1) / 2 + 1)
}
rasterVis::levelplot(result_matrix, att = n,
                     # col.regions = rainbow(odd_number_larger),
                     colorkey = list(at = regions,
                                     # col=rainbow(odd_number_larger)[keep_these_colors],
                                     labels = list(at = middle_of_regions, labels = as.character(1:n))),
                     xlab = 'Cells',
                     ylab = 'Target genes',
                     main='Generated data')

#-------
#----- Construct heatmap for our biclust

plots <- vector(mode = "list", length = length(penalization_lambdas))
for(i_res in 1:length(penalization_lambdas)){
  res <- BICLUST_RESULTS[[i_res]]
  cell_cluster_allocation <- res$cell_cluster_allocation
  target_gene_allocation <- res$scregclust_final_result_module

  result_matrix <- matrix(0, nrow = n_total_cells, ncol = n_target_genes)

  for (i in 1:n_total_cells) {
    cluster <- cell_cluster_allocation[i]
    gene_allocation <- target_gene_allocation[[cluster]][1:n_target_genes]
    gene_allocation[gene_allocation==-1] <- 0
    # Create unique numbers for each pair
    result_matrix[i, ] <- paste0(cluster, gene_allocation)
  }

  # Convert the result to numeric matrix
  result_matrix <- matrix(as.numeric(result_matrix), nrow = n_total_cells, ncol = n_target_genes)
  result_matrix <- matrix(as.integer(as.factor(result_matrix)), nrow=nrow(result_matrix), ncol=ncol(result_matrix))

  # calc_hamming(unique(result_matrix, MARGIN = 1))

  RI <- round(aricode::RI(as.vector(result_matrix), correct_clustering), 2)
  print(RI)

  n <- length(unique(as.vector(result_matrix)))
  regions <- seq(1, n, length.out = n + 1)
  middle_of_regions <- (regions[-1] + regions[-length(regions)]) / 2
  odd_number_larger <- ifelse(n %% 2 == 0, n + 1, n)
  if(n %% 2 == 1){
    keep_these_colors = 1:n
  }else{
    keep_these_colors <- setdiff(1:odd_number_larger, (odd_number_larger + 1) / 2 + 1)
  }
  plots[[i_res]] <- rasterVis::levelplot(result_matrix,
                                         att = n,
                                         # col.regions = rainbow(odd_number_larger),
                                         colorkey = list(at = regions,
                                                         # col=rainbow(odd_number_larger)[keep_these_colors],
                                                         labels = list(at = middle_of_regions, labels = as.character(1:n))),
                                         xlab = 'Cells',
                                         ylab = 'Target genes',
                                         main=paste0('biclust_scregclust, lambda:', penalization_lambdas[i_res], ", RI:", RI))
}

cowplot::plot_grid(plotlist = plots,  align = 'vh', axis = 'tblr')
