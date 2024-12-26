#!/usr/bin/Rscript
rm(list = ls())

library(here)  # To work with paths
library(patchwork)
sink()

# options(warn=2)  # To convert warning messages into error messages which display row of error. For debugging.

# Get absolute path where script is located, by using relative paths.
demo_path <- here::here("demo")
R_path <- here::here("R")
output_path <- demo_path
path_data <- here::here('data')

redo_flag = T

source(file.path(R_path, "generate_dummy_data_for_cell_clustering.R"))
source(file.path(R_path, "bisc.R"))


#############################################
############ data for dev ###################
#############################################

# Set seed for example
set.seed(1234)

# Set variables ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
n_cell_clusters = 2
n_target_gene_clusters = c(4, 2)  # Number of target gene clusters in each cell cluster
n_target_genes = 10
n_regulator_genes = 6
n_cells = c(100, 100)
regulator_means = c(0, 0) # For generating dummy data, regulator mean in each cell cluster
coefficient_means <- list(c(1,3,5,7),              c(10, 20))  # For generating dummy data, coefficient means in each cell cluster
coefficient_sds <- list(c(0.01, 0.01, 0.01, 0.01), c(0.01, 0.01 ))
disturbed_fraction = 0.1  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
plot_stuff = TRUE
plot_suffix = "Simple"
testing_penalization_data_gen = c(0.1, 0.5)


# n_cell_clusters <- 2
# n_target_gene_clusters <- c(4, 2)  # Number of target gene clusters in each cell cluster
# n_target_genes <- 100
# n_regulator_genes <- 10
# n_cells <- c(100, 100)
# regulator_means <- c(0, 0) # For generating dummy data, regulator mean in each cell cluster
# coefficient_means <- list(c(1, 3, 5, 7), c(10, 20))  # For generating dummy data, coefficient means in each cell cluster
# coefficient_sds <- list(c(0.01, 0.01, 0.01, 0.01), c(0.01, 0.01))
# disturbed_fraction <- 0.1  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
# plot_stuff <- FALSE
# plot_suffix <- "Simple"
# testing_penalization_data_gen <- c(0.1, 0.5) #same length as n_cells


if (
      !file.exists(file.path(path_data, "env_sim_simple_data_biclust_sc.rds"))  |  redo_flag
    ) {

  generated_data <- generate_dummy_data_for_cell_clustering(
    n_cell_clusters = 2,
    n_target_gene_clusters = c(4, 2),  # Number of target gene clusters in each cell cluster
    n_target_genes = n_target_genes,          #from vignette
    n_regulator_genes = n_regulator_genes,    # from vignette
    n_cells = n_cells,
    regulator_means = regulator_means,# For generating dummy data, regulator mean in each cell cluster
    coefficient_means <- coefficient_means,  # For generating dummy data, coefficient means in each cell cluster
    coefficient_sds <- coefficient_sds,
    disturbed_fraction = 0.1,  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
    plot_stuff = FALSE,
    plot_suffix = "Simple",
    testing_penalization = testing_penalization_data_gen,
    generate_counts = F,
    trivial_regulator_networks = T,
    pearson_regulators          = T
  )

  saveRDS(generated_data, file.path(path_data, "env_sim_simple_data_biclust_sc.rds"))

} else{

  generated_data <- readRDS(file.path(path_data, "env_sim_simple_data_biclust_sc.rds"))

}

  str(generated_data)

 # generated_data$true_target_gene_allocation

# Because "dat <- cbind(Z_t, Z_r)" in generate_dummy_data_for_cell_clustering
ind_targetgenes <- which(c(rep(1, n_target_genes), rep(0, n_regulator_genes)) == 1)
ind_reggenes <- which(c(rep(0, n_target_genes), rep(1, n_regulator_genes)) == 1)


disturbed_initial_cell_clust <- factor(generated_data$disturbed_initial_cell_clust)

biclust_input_data <- generated_data$dat
colnames(biclust_input_data) <- c(paste0("t", 1:n_target_genes), paste0("r", 1:n_regulator_genes))
biclust_input_data <- tibble::as_tibble(biclust_input_data)

# # These needs to be strings for discrete labels in pca plot
# data_for_plotting <- tibble::as_tibble(true_cell_cluster_allocation = generated_data$true_cell_clust,
#                                        biclust_input_data)
# pca_res <- prcomp(biclust_input_data, scale. = TRUE)}
# p <- ggplot2::autoplot(pca_res, data = data_for_plotting, colour = 'true_cell_cluster_allocation')

# Set up some variables
n_cell_clusters <- length(unique(disturbed_initial_cell_clust))
n_target_genes <- length(ind_targetgenes)
n_regulator_genes <- length(ind_reggenes)


#############################################
############ end data for dev ###############
#############################################

###########################################
####initialise variables for dev ##########
##########################################

max_iter <- 100

# output_path <- modded_output_path

i_cell_cluster <- 1

use_weights <- TRUE
use_complex_cluster_allocation <- FALSE

demo_path <- here::here("demo")
output_path <- demo_path
# initial_clustering[1:14900] <- rep(1, 14900)
# initial_clustering[14901:15000] <- rep(2, 100)
# print(length(initial_clustering))
# stop("hej")


###########################################
###END initialise variables for dev #######
###########################################
# Split data into train/test
cell_data_split <- sample(c(1, 2), nrow(biclust_input_data), prob = c(0.5, 0.5), replace = T)
train_indices <- which(cell_data_split == 1)
test_indices <- which(cell_data_split == 2)

biclust_input_data_train <- biclust_input_data[train_indices,]
biclust_input_data_test <- biclust_input_data[test_indices,]

# Setup variables that will be used throughout the script
# We assume target genes start with t then a number. And likewise for regulator genes.
n_total_cells_train <- length(train_indices)
n_total_cells_test <- length(test_indices)
n_total_cells <- nrow(biclust_input_data)
cell_id <- 1:n_total_cells
cell_id_train <- cell_id[train_indices]
cell_id_test <- cell_id[test_indices]

initial_clustering <- disturbed_initial_cell_clust

initial_clustering_train <- initial_clustering[train_indices]

initial_clustering_test <- initial_clustering[test_indices]

# Set up some variables

# true_cell_cluster_allication <- factor(generated_data$true_cell_clust)
# true_cell_cluster_allication_train <- true_cell_cluster_allication[train_indices]
# true_cell_cluster_allication_train <- true_cell_cluster_allication[train_indices]


penalization_lambdas <- c( 0.1, 0.2, 0.5) # c( 0.00001, 0.1, 0.2, 0.5)
BICLUST_RESULTS_train <- vector(mode = "list", length = length(penalization_lambdas))

max_iter <- 20

if (
  !file.exists(file.path(path_data, "env_sim_simple_nogarb_res_biclust_sc.rds"))  | redo_flag
) {

  for (i_penalization_lambda in seq_along(penalization_lambdas)) {
    print("", quote = FALSE)
    print(paste("Running biclust for penalization_lambda", penalization_lambdas[i_penalization_lambda]), quote = FALSE)
    BICLUST_RESULTS_train[[i_penalization_lambda]] <- bisc(
      dat = biclust_input_data_train,
      cell_id = cell_id_train,
      true_cell_cluster_allocation = factor(generated_data$true_cell_clust[train_indices]),
      max_iter = max_iter,
      n_target_gene_clusters,
      initial_clustering_train,
      n_cell_clusters,
      ind_targetgenes,
      ind_reggenes,
      output_path,
      penalization_lambda = penalization_lambdas[i_penalization_lambda],
      calculate_optimization_target = TRUE,
      calculate_silhoutte = FALSE,
      calculate_davies_bouldin_index = TRUE,
      plot_suffix = "Simple_cluster_all",
      always_use_flat_prior = FALSE,
      use_garbage_cluster_targets  = F
    )
  }

  saveRDS(BICLUST_RESULTS_train, file.path(path_data, "env_sim_simple_nogarb_res_biclust_sc.rds"))

} else{

  BICLUST_RESULTS_train <- readRDS(file.path(path_data, "env_sim_simple_nogarb_res_biclust_sc.rds"))

}

str(BICLUST_RESULTS_train)


print("", quote = FALSE)
print("", quote = FALSE)
for (i_penalization_lambda in seq_along(penalization_lambdas)) {
  if (is.na(BICLUST_RESULTS_train[i_penalization_lambda])) {
    print(paste("penalization_lambda", penalization_lambdas[i_penalization_lambda], "is NA"), quote = FALSE)
  } else if (is.null(BICLUST_RESULTS_train[i_penalization_lambda])) {
    print(paste("penalization_lambda", penalization_lambdas[i_penalization_lambda], "is NULL"), quote = FALSE)
  }else {
    print(paste("penalization_lambda", penalization_lambdas[i_penalization_lambda], "is ok with rand index", BICLUST_RESULTS_train[[i_penalization_lambda]]$rand_index), quote = FALSE)
  }
}


################################################################
############### what do we plot ??? ############################
################################################################


RIs <- sapply(BICLUST_RESULTS_train, function(x) as.numeric(x$rand_index))





png(file.path(output_path, paste0("bisc_lambda_RI_simple_ex",".png")),
    width = 1024, height = 480, units = "px")
plot(
  penalization_lambdas,
  RIs,
  type = 'l'
)
dev.off()



############################################################################################
############ run the penalisation parameter by best BIC, on the test dataset#################

bic_values <- sapply(BICLUST_RESULTS_train, function(x) as.numeric(x$BIC[[length(x$BIC)]]))

print(bic_values)

best_lambda <- penalization_lambdas[which.min(bic_values)] # pick best
BICLUST_RESULTS_test <- vector(mode = "list", length = length(penalization_lambdas))

max_iter <- 20

if (
  !file.exists(file.path(path_data, "env_sim_simple_nogarb_res_biclust_sc_test.rds"))  | redo_flag
) {

  for (i_penalization_lambda in best_lambda) {
    print("", quote = FALSE)
    BICLUST_RESULTS_test<- bisc(
      dat = biclust_input_data_test,
      cell_id = cell_id_test,
      true_cell_cluster_allocation = factor(generated_data$true_cell_clust[test_indices]),
      max_iter = max_iter,
      n_target_gene_clusters,
      initial_clustering_test,
      n_cell_clusters,
      ind_targetgenes,
      ind_reggenes,
      output_path,
      penalization_lambda = best_lambda,
      calculate_optimization_target = TRUE,
      calculate_silhoutte = FALSE,
      calculate_davies_bouldin_index = TRUE,
      plot_suffix = "Simple_cluster_all",
      always_use_flat_prior = FALSE,
      use_garbage_cluster_targets  = F
    )
  }

  saveRDS(BICLUST_RESULTS_test, file.path(path_data, "env_sim_simple_nogarb_res_biclust_sc_test.rds"))

} else{

  BICLUST_RESULTS_test <- readRDS(file.path(path_data, "env_sim_simple_nogarb_res_biclust_sc_test.rds"))

}

####################################
#### what stats/ results do we want?
#######################################



##############################
# run again many times for stats
###########################


iterations <- 20

BICLUST_RESULTS_iterated <- vector(mode = "list", length = iterations)


if (
  !file.exists(file.path(path_data, "env_sim_simple_nogarb_res_biclust_sc_iter.rds"))  | redo_flag
) {

  for (i_penalization_lambda in 1:iterations) {
    print(i_penalization_lambda, quote = FALSE)
    BICLUST_RESULTS_iterated[[i_penalization_lambda]] <- bisc(
      dat = biclust_input_data_test,
      cell_id = cell_id_test,
      true_cell_cluster_allocation = factor(generated_data$true_cell_clust[test_indices]),
      max_iter = max_iter,
      n_target_gene_clusters,
      initial_clustering_test,
      n_cell_clusters,
      ind_targetgenes,
      ind_reggenes,
      output_path,
      penalization_lambda = best_lambda,
      calculate_optimization_target = TRUE,
      calculate_silhoutte = FALSE,
      calculate_davies_bouldin_index = TRUE,
      plot_suffix = "Simple_cluster_all",
      always_use_flat_prior = FALSE,
      use_garbage_cluster_targets  = F
    )
  }

  saveRDS(BICLUST_RESULTS_iterated, file.path(path_data, "env_sim_simple_nogarb_res_biclust_sc_iter.rds"))

} else{

  BICLUST_RESULTS_iterated <- readRDS(file.path(path_data, "env_sim_simple_nogarb_res_biclust_sc_iter.rds"))

}

str(BICLUST_RESULTS_iterated[[1]])
BICLUST_RESULTS_iterated[[2]]$cell_cluster_allocation %>% table()
factor(generated_data$true_cell_clust[test_indices]) %>% table()
# RI(unlist(BICLUST_RESULTS_iterated[[2]]$cell_cluster_allocation),factor(generated_data$true_cell_clust[test_indices]) )

all(is.na(BICLUST_RESULTS_iterated[[2]]))

BICLUST_RESULTS_iterated[[2]]$rand_index
#make some plots of results
RIs <-  sapply(seq_along(BICLUST_RESULTS_iterated),
               function(x) {
                 print(x)
                 if(all(is.na(BICLUST_RESULTS_iterated[[x]]))){
                   NA
                 }else{
                   as.numeric(BICLUST_RESULTS_iterated[[x]]$rand_index)
                 }
               }
               )


png(file.path(output_path, paste0("bisc_iterated_RI_histogram_simple_ex",".png")),
    width = 1024, height = 480, units = "px")
hist( RIs,
      breaks = 10,
      main = paste0("Historam of RIs, \n",iterations, " iterations\n", "Mean RI: ", round(mean(RIs), 3)))

dev.off()


#start by using biclust to get som cell clusters only, and use that one for plots, update later.
#
# # compare with other biclust variant
# library(biclust)
#
# # ?biclust
#
# # test <- matrix(rbinom(400, 50, 0.4), 20, 20)
# # res1 <- biclust(as.matrix(biclust_input_data), method=BCCC(), delta=1.5,  alpha=1, number=10)
# # res2 <- biclust(as.matrix(biclust_input_data), method=BCPlaid(),  cluster="b")
#
#
# comparator_iterated <- vector(mode = "list", length = iterations)
#
# if (
#   !file.exists(file.path(path_data, "env_sim_simple_comparator_iter.rds"))  | redo_flag
# ) {
#
#   for (iter in 1:iterations) {
#     print(iter)
#     comparator_iterated[[iter]] <- biclust(as.matrix(biclust_input_data), method=BCPlaid(),  cluster="r")
#
#   }
#
#   saveRDS(comparator_iterated, file.path(path_data, "env_sim_simple_comparator_iter.rds"))
#
#   } else {
#
#     comparator_iterated <- readRDS(file.path(path_data, "env_sim_simple_comparator_iter.rds"))
#
# }
#
# str(comparator_iterated[[1]])
#
# clusters_comparator <-  lapply(comparator_iterated, function(x) x@RowxNumber )
#
#
#
# sapply(1:nrow(clusters_comparator), FUN = function(i) which(clusters_comparator@RowxNumber[i,]))
#
#
# str(res2)
# res2@RowxNumber
# sort(unique(unlist(sapply(1:nrow(biclust_input_data), FUN = function(i) which(res2@RowxNumber[i,])))))
#
#
# res1@RowxNumber
#
# image(res1@RowxNumber)
# ?biclust::heatmapBC


# Celda variant?
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#
# BiocManager::install("celda", force = TRUE)
#
# browseVignettes("celda")

library("celda")

test <- SingleCellExperiment(  assays = list ( counts = t(as.matrix(generated_data$counts
                                                                    ) + 1
                                                          )
                                              )
                               )
str(test)
rowData(test)
colData(test)

test <- selectFeatures(test)
sce <- celda_CG(x = test , # cells should be columns here
                K = 2, L = 4, verbose = TRUE, nchains = 1)

table(celdaClusters(sce))
table(celdaModules(sce))
sce <- celdaUmap(sce)

plotDimReduceCluster(x = sce, reducedDimName = "celda_UMAP")
celdaProbabilityMap(sce)
plot(celdaHeatmap(sce = sce, nfeatures = 10))

moduleHeatmap(sce, featureModule = c(1,2), topCells = 100)


###############################################################################
################## try running the counts through sctransform #################
################## nd see if biclust does anything useful######################

library(Seurat)
library(sctransform)
simdata <- CreateSeuratObject(counts = t(generated_data$counts))

# store mitochondrial percentage in object meta data
simdata <- PercentageFeatureSet(simdata, pattern = "^MT-", col.name = "percent.mt")

# run sctransform
simdata <- SCTransform(simdata, vars.to.regress = "percent.mt", verbose = FALSE)


simdata <- RunPCA(simdata, verbose = FALSE)
simdata <- RunUMAP(simdata, dims = 1:30, verbose = FALSE)

simdata <- FindNeighbors(simdata, dims = 1:30, verbose = FALSE)
simdata <- FindClusters(simdata, verbose = FALSE)
DimPlot(simdata, label = TRUE)

simdata@assays$SCT@SCTModel.list$counts@feature.attributes$theta

z <- GetAssayData(simdata, layer = "scale.data")
dim(z)


# print(paste("rand_index for result vs true cluster:", biclust_result$rand_index), quote = FALSE)
# print(paste("Number of iterations:", biclust_result$n_iterations), quote = FALSE)
# print(paste("Silhoutte of first disturbed cluster likelihood (aka how complex was the first likelihood):", biclust_result$db), quote = FALSE)
# print(paste("BIC_all:", biclust_result$BIC), quote = FALSE)
# print(paste("target_genes_residual_var:"), quote = FALSE)
# print(biclust_result$taget_genes_residual_var, quote = FALSE)
