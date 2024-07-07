#!/usr/bin/Rscript
rm(list = ls())

library(here)  # To work with paths
library(patchwork)
library(scregclust)
library("EnsDb.Hsapiens.v79")  # BiocManager::install("EnsDb.Hsapiens.v79")
sink()
gc()  # Force clean memory

# options(warn=2)  # To convert warning messages into error messages which display row of error. For debugging.

# Get absolute path where script is located, by using relative paths.
demo_path <- here::here("demo")
R_path <- here::here("R")
output_path <- demo_path


#############################################
############ data for dev ###################
#############################################

# Set seed for example
set.seed(214)

# Load data
# Folders and filenames
path_data <- here::here('data')
path_medium <- file.path(path_data, "medium.rds")
path_medium_cleaned <- file.path(path_data, "medium_cleaned_findtargetcluster.rds")
path_medium_cleaned_sctransform <- file.path(path_data, "medium_cleaned_sctransform_findtargetcluster.rds")

path_cell_types <- file.path(path_data, "cell_types_findtargetcluster.rds")
path_env_data <- file.path(path_data, "env_data_findtargetcluster_sctransform.RData")
path_general_env_data <- file.path(path_data, "env_data_general_findtargetcluster_sctransform.RData")
path_AC <- file.path(path_data, "env_data_AC_sctransform.RData")
path_MES <- file.path(path_data, "env_data_MES_sctransform.RData")
path_NPC <- file.path(path_data, "env_data_NPC_sctransform.RData")
path_OPC <- file.path(path_data, "env_data_OPC_sctransform.RData")

path_env_data_s <- file.path(path_data, "env_data_findtargetcluster_sctransform_s.RData")
path_general_env_data_s <- file.path(path_data, "env_data_general_findtargetcluster_sctransform_s.RData")
path_AC_s <- file.path(path_data, "env_data_AC_sctransform_s.RData")
path_MES_s <- file.path(path_data, "env_data_MES_sctransform_s.RData")
path_NPC_s <- file.path(path_data, "env_data_NPC_sctransform_s.RData")
path_OPC_s <- file.path(path_data, "env_data_OPC_sctransform_s.RData")
# Data from https://cellxgene.cziscience.com/collections/999f2a15-3d7e-440b-96ae-2c806799c08c

path_env_data_tsne3 <- file.path(path_data, "env_data_findtargetcluster_sctransform_tsne3.RData")

if (file.exists(path_env_data)) {
  load(path_env_data_s)
}else {
  if (file.exists(path_medium_cleaned) && file.exists(path_cell_types)) {
    # Read data
    d <- readRDS(path_medium_cleaned)
    cell_types <- readRDS(path_cell_types)
  }else {
    # Read data
    d <- readRDS(path_medium)

    # Keep the correct type of cells
    logical_keep_neftel_cells <- as.character(d@meta.data$author)=="Neftel2019"
    cell_types <- d@meta.data$annotation_level_3
    keep_cells <- cell_types == 'AC-like' |
      cell_types == 'MES-like' |
      cell_types == 'NPC-like' |
      cell_types == 'OPC-like'
    keep_cells <- logical_keep_neftel_cells & keep_cells


    cell_types <- factor(cell_types[keep_cells])
    saveRDS(cell_types, path_cell_types)
    d <- d[, keep_cells]

    # Remove genes which are 'constant', aka no variance, just one number.
    # scregclust doesn't work on those.
    d_temp <- d[["RNA"]]$data
    non_constant_ind <- which(apply(d_temp, 1, sd) != 0)
    d <- d[non_constant_ind,]
    rm(non_constant_ind)
    rm(d_temp)
    saveRDS(d, path_medium_cleaned)
  }

  # Find out which are regulator genes
  # Make it into a matrix and rename row names from
  # Ensembl gene IDs (like "ENSG00000269696") to standard gene names (like "A1BG"),
  fake_matrix <- matrix(0, nrow = nrow(d), ncol = 1)  # This is fast
  # dm <- tibble::as_tibble(as.matrix(GetAssayData(d, assay = "RNA", slot = "counts")))  # Using the data is slow
  rownames(fake_matrix) <- as.vector(d@assays$RNA@meta.features$feature_name)
  out <- scregclust::scregclust_format(fake_matrix)  # Needs to be a matrix to use this
  is_regulator <- out[['is_regulator']]
  rm(fake_matrix, out)

  # Put some stuff in order gene/row wise.
  # This puts the regulator genes at the end
  d <- d[order(is_regulator),]
  is_regulator <- is_regulator[order(is_regulator)]
  n_total_regulators <- sum(is_regulator)
  n_total_target_genes <- nrow(d) - n_total_regulators

  print(object.size(d), units = 'MB', standard = 'SI')

  # Sub sample genes
  n_target_genes <- n_total_target_genes  # 25518
  n_regulator_genes <- n_total_regulators  # 1562
  if (n_target_genes < n_total_target_genes || n_regulator_genes < n_total_regulators) {
    old_ind_regulators <- which(as.logical(is_regulator))
    new_ind_targetgenes <- sample(1:(min(old_ind_regulators) - 1), n_target_genes)
    new_ind_regulators <- sample(old_ind_regulators, n_regulator_genes)
    new_ind_genes <- sort(c(new_ind_targetgenes, new_ind_regulators))
    keep_genes <- replace(rep(FALSE, length(is_regulator)), new_ind_genes, TRUE)
    #  Apply sub sampling
    d <- d[keep_genes,]
    is_regulator <- (c(rep(0, n_target_genes), rep(1, n_regulator_genes)) == 1)
    rm(keep_genes, new_ind_genes, new_ind_regulators, new_ind_targetgenes, old_ind_regulators)
    n_target_genes <- length(is_regulator) - sum(is_regulator == 1)
    n_regulator_genes <- sum(is_regulator == 1)
  }
  rm(n_total_regulators, n_total_target_genes)

  # Put stuff in order cell/column wise
  cell_order <- order(cell_types)
  cell_types <- cell_types[cell_order]
  d <- d[, cell_order]
  rm(cell_order)

  # Now when we have put the cells in order we just need to count the cells
  # in each cell cluster. Then it's easy to make a vector with the true_cluster_allocation
  # further down in code
  n_cells_cell_cluster_1 <- sum(cell_types == 'AC-like')
  n_cells_cell_cluster_2 <- sum(cell_types == 'MES-like')
  n_cells_cell_cluster_3 <- sum(cell_types == 'NPC-like')
  n_cells_cell_cluster_4 <- sum(cell_types == 'OPC-like')

  print(paste("Cells in clusters:"), quote = FALSE)
  print(paste("AC-like:", n_cells_cell_cluster_1), quote = FALSE)
  print(paste("MES-like:", n_cells_cell_cluster_2), quote = FALSE)
  print(paste("NPC-like:", n_cells_cell_cluster_3), quote = FALSE)
  print(paste("OPC-like:", n_cells_cell_cluster_4), quote = FALSE)
  print(paste("Number of regulator genes are", sum(is_regulator)), quote = FALSE)
  print("Scregclust wants more cells than regulator genes x 2 for each cell cluster. Otherwise it doesn't work.", quote = FALSE)
  rm(cell_types)

  # Set variables ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  n_cell_clusters <- 4
  # We assume cells are ordered in the order of cell clusters. So the first x columns are cell cluster 1, etc.
  n_cells_in_each_cluster <- c(n_cells_cell_cluster_1, n_cells_cell_cluster_2, n_cells_cell_cluster_3, n_cells_cell_cluster_4)
  true_cluster_allocation <- rep(1:n_cell_clusters, times = n_cells_in_each_cluster)  # TODO: do this
  rm(n_cells_in_each_cluster, n_cells_cell_cluster_1, n_cells_cell_cluster_2, n_cells_cell_cluster_3, n_cells_cell_cluster_4)

  ind_targetgenes <- which(c(rep(1, n_target_genes), rep(0, n_regulator_genes)) == 1)
  ind_reggenes <- which(c(rep(0, n_target_genes), rep(1, n_regulator_genes)) == 1)

  # Convert from seurat object into matrix that we can work with
  # List of object layers
  # SeuratObject::Layers(d)
  #-------------------------------------------------------------------
  d <- Seurat::SCTransform(d,variable.features.n = 4000)
  # d <- SCTransform(d, variable.features.n = 27080)
  # saveRDS(d, path_medium_cleaned_sctransform)
  # Find out which are regulator genes
  # Make it into a matrix and rename row names from
  # Ensembl gene IDs (like "ENSG00000269696") to standard gene names (like "A1BG"),
  ensemblgenes <- d@assays$SCT@var.features
  geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensemblgenes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
  select_these <- geneIDs1[,2] # Since we probably lost genes in the translation
  new_names <- geneIDs1[,1]


  d <- d[select_these,]
  fake_matrix <- matrix(0, nrow = nrow(d), ncol = 1)  # This is fast
  # dm <- tibble::as_tibble(as.matrix(GetAssayData(d, assay = "RNA", slot = "counts")))  # Using the data is slow
  rownames(fake_matrix) <- new_names  # as.vector(d@assays$SCT@var.features)
  out <- scregclust::scregclust_format(fake_matrix)  # Needs to be a matrix to use this
  is_regulator <- out[['is_regulator']]
  rm(fake_matrix, out)



  #------------------------------------------------------------------
  biclust_input_data <- as.matrix(Seurat::GetAssayData(d, assay = "SCT", layer = "scale.data"))
  org_gene_names <- rownames(biclust_input_data)
  org_cell_names <- colnames(biclust_input_data)
  rm(d)



  biclust_input_data <- t(tibble::as_tibble(biclust_input_data))  # TODO: WHYYYY. It takes more memory.

  n_regulator_genes <- sum(is_regulator)
  n_target_genes <- ncol(biclust_input_data) - n_regulator_genes
  new_gene_names <- c(paste0("t", 1:n_target_genes), paste0("r", 1:n_regulator_genes))
  ind_targetgenes <- which(c(rep(1, n_target_genes), rep(0, n_regulator_genes)) == 1)
  ind_reggenes <- which(c(rep(0, n_target_genes), rep(1, n_regulator_genes)) == 1)

  colnames(biclust_input_data) <- new_gene_names
  rm(new_gene_names)

  gc()  # Force clean memory

  save(biclust_input_data,
       org_gene_names,
       org_cell_names,
       n_cell_clusters,
       n_regulator_genes,
       n_target_genes,
       ind_reggenes,
       ind_targetgenes,
       true_cluster_allocation,
       file = path_env_data_s)

   save(n_cell_clusters,
       org_gene_names,
       org_cell_names,
       n_regulator_genes,
       n_target_genes,
       ind_reggenes,
       ind_targetgenes,
       true_cluster_allocation,
       file = path_general_env_data_s)
}

if (!file.exists(path_AC_s) ||
  !file.exists(path_MES_s) ||
  !file.exists(path_NPC_s) ||
  !file.exists(path_OPC_s)) {
  print(object.size(biclust_input_data), units = 'MB', standard = 'SI')


  # Remove all genes with less non zero numbers than 200
  biclust_input_data <- biclust_input_data[, apply(biclust_input_data, MARGIN = 2, function(x) sum(x != 0)) > 200]
  n_target_genes <- length(grep("^t", colnames(biclust_input_data), value = TRUE))
  n_regulator_genes <- length(grep("^r", colnames(biclust_input_data), value = TRUE))
  ind_targetgenes <- which(c(rep(1, n_target_genes), rep(0, n_regulator_genes)) == 1)
  ind_reggenes <- which(c(rep(0, n_target_genes), rep(1, n_regulator_genes)) == 1)
  colnames(biclust_input_data) <- c(paste0("t", 1:n_target_genes), paste0("r", 1:n_regulator_genes))
  print(object.size(biclust_input_data), units = 'MB', standard = 'SI')

  load(path_general_env_data_s)

  # ------------------
  # Save one variable per cell cluster to save on memory use
  cell_cluster_AC <- t(biclust_input_data[true_cluster_allocation == 1,])
  print(dim(cell_cluster_AC))
  print(object.size(cell_cluster_AC), units = 'MB', standard = 'SI')
  save(cell_cluster_AC, file = path_AC_s)
  rm(cell_cluster_AC)
  gc()

  cell_cluster_MES <- t(biclust_input_data[true_cluster_allocation == 2,])
  print(dim(cell_cluster_MES))
  print(object.size(cell_cluster_MES), units = 'MB', standard = 'SI')
  save(cell_cluster_MES, file = path_MES_s)
  rm(cell_cluster_MES)
  gc()

  cell_cluster_NPC <- t(biclust_input_data[true_cluster_allocation == 3,])
  print(dim(cell_cluster_NPC))
  print(object.size(cell_cluster_NPC), units = 'MB', standard = 'SI')
  save(cell_cluster_NPC, file = path_NPC_s)
  rm(cell_cluster_NPC)
  gc()

  cell_cluster_OPC <- t(biclust_input_data[true_cluster_allocation == 4,])
  print(dim(cell_cluster_OPC))
  print(object.size(cell_cluster_OPC), units = 'MB', standard = 'SI')
  save(cell_cluster_OPC, file = path_OPC_s)
  rm(cell_cluster_OPC)
  gc()
}


# TODO: Start here
# ----------------------------
load(path_general_env_data_s)
min_number_of_clusters <- 2
max_number_of_clusters <- 8
penalization_lambda <- 2
rm(biclust_input_data)
gc()
target_gene_cluster_vector <- seq(min_number_of_clusters, max_number_of_clusters)

ensemble_to_normal <- function(genes = org_gene_names){
  ensemblgenes <- genes
  geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensemblgenes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
  select_these <- geneIDs1[,2] # Since we probably lost genes in the translation
  normal_gene_names <- geneIDs1[,1]
  return(normal_gene_names)
}

give_me_regulators <- function(normal_gene_names) {
  fake_matrix <- matrix(0, nrow = length(normal_gene_names), ncol = 1)  # This is fast
  rownames(fake_matrix) <- normal_gene_names  # as.vector(d@assays$SCT@var.features)
  out <- scregclust::scregclust_format(fake_matrix)  # Needs to be a matrix to use this
  is_regulator <- out[['is_regulator']]
  rm(fake_matrix, out)
  return(is_regulator)
}

for (i_cell_cluster in seq(n_cell_clusters)) {
  gc()
  print(i_cell_cluster)
  if(i_cell_cluster==1){
    load(path_AC_s)
    current_cell_cluster <- cell_cluster_AC
    rm(cell_cluster_AC)
  }else if(i_cell_cluster==2){
    load(path_MES_s)
    current_cell_cluster <- cell_cluster_MES
    rm(cell_cluster_MES)
  }else if(i_cell_cluster==3){
    load(path_NPC_s)
    current_cell_cluster <- cell_cluster_NPC
    rm(cell_cluster_NPC)
  }else if(i_cell_cluster==4){
    load(path_OPC_s)
    current_cell_cluster <- cell_cluster_OPC
    rm(cell_cluster_OPC)
  }
  sink()
  print("Loading compelte")
  normal_gene_names <- ensemble_to_normal(genes=org_gene_names)
  is_regulator <- give_me_regulators(normal_gene_names)
  print("Running reg")
  sink()
  tsne_3_reg <- Rtsne::Rtsne(t(current_cell_cluster[is_regulator, ]), dims=3, num_threads = 30, verbose=TRUE, partial_pca=TRUE, max_iter=2000, initial_dims=200)
  tsne_3_target <- Rtsne::Rtsne(t(current_cell_cluster[!is_regulator, ]), dims=3, num_threads = 30, verbose=TRUE, partial_pca=TRUE, max_iter=2000, initial_dims=200)
  print("Running target")
  sink()
  full_3t_3r <- t(cbind(tsne_3_reg$Y, tsne_3_target$Y))
    ind_targetgenes <- 1:3
  ind_reggenes <- 4:6
  n_target_genes <- 3
  n_regulator_genes <- 3
  rownames(full_3t_3r) <- c(paste0("t", 1:n_target_genes), paste0("r", 1:n_regulator_genes))

   if(i_cell_cluster==1){
    cell_cluster_AC <- full_3t_3r
  }else if(i_cell_cluster==2){
    cell_cluster_MES <- full_3t_3r
  }else if(i_cell_cluster==3){
    cell_cluster_NPC <- full_3t_3r
  }else if(i_cell_cluster==4){
    cell_cluster_OPC <- full_3t_3r
  }
  rm(full_3t_3r, tsne_3_reg, tsne_3_target)
}

  save(cell_cluster_AC,
       cell_cluster_MES,
       cell_cluster_NPC,
       cell_cluster_OPC,
       org_gene_names,
       org_cell_names,
       n_cell_clusters,
       n_regulator_genes,
       n_target_genes,
       ind_reggenes,
       ind_targetgenes,
       true_cluster_allocation,
       file = path_env_data_tsne3)

rm(max_number_of_clusters, min_number_of_clusters, normal_gene_names, org_gene_names, penalization_lambda, target_gene_cluster_vector)

for (i_cell_cluster in seq(n_cell_clusters)) {
  if(i_cell_cluster==1){
    load(path_AC_s)
    current_cell_cluster <- cell_cluster_AC
    rm(cell_cluster_AC)
  }else if(i_cell_cluster==2){
    load(path_MES_s)
    current_cell_cluster <- cell_cluster_MES
    rm(cell_cluster_MES)
  }else if(i_cell_cluster==3){
    load(path_NPC_s)
    current_cell_cluster <- cell_cluster_NPC
    rm(cell_cluster_NPC)
  }else if(i_cell_cluster==4){
    load(path_OPC_s)
    current_cell_cluster <- cell_cluster_OPC
    rm(cell_cluster_OPC)
  }

  # Eemove cells and target genes
  keep_cells_split <- sample(c(1, 2), ncol(current_cell_cluster), prob = c(1.0, 0), replace = T)
  keep_ind_cells <- keep_cells_split == 1

  current_cell_cluster <- current_cell_cluster[, keep_ind_cells]

  ensemblgenes <- org_gene_names
  geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensemblgenes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
  select_these <- geneIDs1[,2] # Since we probably lost genes in the translation
  new_names <- geneIDs1[,1]

  fake_matrix <- matrix(0, nrow = nrow(current_cell_cluster), ncol = 1)  # This is fast
  # dm <- tibble::as_tibble(as.matrix(GetAssayData(d, assay = "RNA", slot = "counts")))  # Using the data is slow
  rownames(fake_matrix) <- new_names  # as.vector(d@assays$SCT@var.features)
  out <- scregclust::scregclust_format(fake_matrix)  # Needs to be a matrix to use this
  is_regulator <- out[['is_regulator']]
  rm(fake_matrix, out)

  n_regulator_genes <- sum(is_regulator)
  n_target_genes <- nrow(current_cell_cluster) - n_regulator_genes
  ind_targetgenes <- which(c(rep(1, n_target_genes), rep(0, n_regulator_genes)) == 1)
  ind_reggenes <- which(c(rep(0, n_target_genes), rep(1, n_regulator_genes)) == 1)

  # Remove target genes
  keep_target_genes_split <- sample(c(1, 2), n_target_genes, prob = c(1.0, 0.0), replace = T)
  keep_ind_target_genes <- keep_target_genes_split == 1
  ind_targetgenes <- ind_targetgenes[keep_ind_target_genes]
  n_target_genes <- length(ind_targetgenes)
  ind_keep_genes <- c(ind_targetgenes, ind_reggenes)

  current_cell_cluster <- current_cell_cluster[ind_keep_genes, ]

  ind_targetgenes <- which(c(rep(1, n_target_genes), rep(0, n_regulator_genes)) == 1)
  ind_reggenes <- which(c(rep(0, n_target_genes), rep(1, n_regulator_genes)) == 1)

  fake_matrix <- matrix(0, nrow = nrow(current_cell_cluster), ncol = 1)  # This is fast
  # dm <- tibble::as_tibble(as.matrix(GetAssayData(d, assay = "RNA", slot = "counts")))  # Using the data is slow
  rownames(fake_matrix) <- new_names[ind_keep_genes]  # as.vector(d@assays$SCT@var.features)
  out <- scregclust::scregclust_format(fake_matrix)  # Needs to be a matrix to use this
  is_regulator <- out[['is_regulator']]
  rm(fake_matrix, out)

  # Remove regulator genes
  keep_regulator_genes_split <- sample(c(1, 2), n_regulator_genes, prob = c(1.0, 0.0), replace = T)
  keep_ind_regulator_genes <- keep_regulator_genes_split == 1
  ind_reggenes <- ind_reggenes[keep_ind_regulator_genes]
  n_regulator_genes <- length(ind_reggenes)
  ind_keep_genes <- c(ind_targetgenes, ind_reggenes)

  current_cell_cluster <- current_cell_cluster[ind_keep_genes, ]

  ind_targetgenes <- which(c(rep(1, n_target_genes), rep(0, n_regulator_genes)) == 1)
  ind_reggenes <- which(c(rep(0, n_target_genes), rep(1, n_regulator_genes)) == 1)

  fake_matrix <- matrix(0, nrow = nrow(current_cell_cluster), ncol = 1)  # This is fast
  # dm <- tibble::as_tibble(as.matrix(GetAssayData(d, assay = "RNA", slot = "counts")))  # Using the data is slow
  rownames(fake_matrix) <- new_names[ind_keep_genes]  # as.vector(d@assays$SCT@var.features)
  out <- scregclust::scregclust_format(fake_matrix)  # Needs to be a matrix to use this
  is_regulator <- out[['is_regulator']]
  rm(fake_matrix, out)


  # Run screg with a bunch of different cluster setings
  results <- vector(mode = "list", length = length(target_gene_cluster_vector))

  for (i_n_target_genes_clusters in seq(length(target_gene_cluster_vector))) {
    n_target_genes_clusters <- target_gene_cluster_vector[i_n_target_genes_clusters]
    print(paste("Cell cluster", i_cell_cluster, "Number of target gene clusters",n_target_genes_clusters), quote=FALSE)

    results[[i_n_target_genes_clusters]] <- scregclust::scregclust(
      expression = current_cell_cluster,  # p rows of genes, n columns of cells
      split_indices = NULL,
      genesymbols = 1:(n_target_genes + n_regulator_genes),  # Gene row numbers
      is_regulator = is_regulator, # inverse_which(indices = ind_reggenes, output_length = n_regulator_genes + n_target_genes),  # Vector indicating which genes are regulators
      n_cl = n_target_genes_clusters,
      penalization = 0.000001,
      noise_threshold = 0.000001,
      verbose = TRUE,
      n_cycles = 40,
      compute_silhouette = TRUE,
      center=TRUE
    )
  }

  saveRDS(results, file.path(path_data, paste0("screg_results_lambda_", penalization_lambda,
                                               "_cell_cluster_", i_cell_cluster,
                                               "_reg_", n_regulator_genes,
                                               "_target_", n_target_genes,
                                               "_cells_", ncol(current_cell_cluster),
                                               "_mincluster_", min_number_of_clusters,
                                                "_maxcluster", max_number_of_clusters,".rds")))
  rm(results, current_cell_cluster)
  gc()
}



scregclust::plot_silhouettes(list_of_fits = results,
                             penalization = penalization_lambda)

scregclust::plot_cluster_count_helper(list_of_fits = results, penalization = penalization_lambda)
penalization_lambda <- c(0.1, 0.2, 0.3, 0.5)
path_data <- "C:\\Users\\Sebastian\\repos\\biclust\\demo"
results1 <- readRDS(file.path(path_data,"scregResultsNeftel2019_lambda_0.1-0.2-0.3-0.5_cellCluster_1_nRegulatorGenes_208_nTargetGenes_1640_nCells_4872_minNCluster_2_maxNCluster15.rds"))
results2 <- readRDS(file.path(path_data,"scregResultsNeftel2019_lambda_0.1-0.2-0.3-0.5_cellCluster_2_nRegulatorGenes_208_nTargetGenes_1640_nCells_2923_minNCluster_2_maxNCluster15.rds"))
results3 <- readRDS(file.path(path_data,"scregResultsNeftel2019_lambda_0.1-0.2-0.3-0.5_cellCluster_3_nRegulatorGenes_208_nTargetGenes_1640_nCells_2810_minNCluster_2_maxNCluster15.rds"))
results4 <- readRDS(file.path(path_data,"scregResultsNeftel2019_lambda_0.1-0.2-0.3-0.5_cellCluster_4_nRegulatorGenes_208_nTargetGenes_1640_nCells_2093_minNCluster_2_maxNCluster15.rds"))
names <- list("AC", "MES", "NPS", "OPC")
all_results <- list(results1, results2, results3, results4)
for(i in seq(4)){
  results <- all_results[i]
  p <- scregclust::plot_cluster_count_helper(list_of_fits = results[[1]], penalization = penalization_lambda)
  print(p)
}


for (i_cell_cluster in seq(n_cell_clusters)) {
  if(i_cell_cluster==1){
    load(path_AC_s)
    current_cell_cluster <- cell_cluster_AC
    print(dim(current_cell_cluster))
    rm(cell_cluster_AC)
  }else if(i_cell_cluster==2){
    load(path_MES_s)
    current_cell_cluster <- cell_cluster_MES
     print(dim(current_cell_cluster))
    rm(cell_cluster_MES)
  }else if(i_cell_cluster==3){
    load(path_NPC_s)
    current_cell_cluster <- cell_cluster_NPC
     print(dim(current_cell_cluster))
    rm(cell_cluster_NPC)
  }else if(i_cell_cluster==4){
    load(path_OPC_s)
    current_cell_cluster <- cell_cluster_OPC
     print(dim(current_cell_cluster))
    rm(cell_cluster_OPC)
  }
}
