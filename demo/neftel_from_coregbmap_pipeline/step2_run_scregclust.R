#!/usr/bin/Rscript
rm(list = ls())

library(here)  # To work with paths
# library(patchwork)
library(scregclust)
# To be able to run library("EnsDb.Hsapiens.v79") you need to at least:
# BiocManager::install("GenomeInfoDb")
# BiocManager::install("SparseArray")
# BiocManager::install("EnsDb.Hsapiens.v79")
library("EnsDb.Hsapiens.v79")
sink()  # Because some of our scripts redirects output from scregclust to force it to be quite. This restores output.
gc()  # Force clean memory

# options(warn=2)  # To convert warning messages into error messages which display row of error. For debugging.

# Get absolute path where script is located, by using relative paths.
demo_path <- here::here("demo")
R_path <- here::here("R")
output_path <- demo_path

path_data <- here::here('data')

# Data from https://cellxgene.cziscience.com/collections/999f2a15-3d7e-440b-96ae-2c806799c08c
path_coreGBmap <- file.path(path_data, "medium.rds")
path_env_data_neftel2019 <- file.path(path_data, "env_data_findtargetcluster_sctransform_neftel2019.RData")
path_general_env_data_neftel2019 <- file.path(path_data, "env_data_general_findtargetcluster_sctransform_neftel2019.RData")
path_AC_neftel2019 <- file.path(path_data, "env_data_AC_sctransform_neftel2019.RData")
path_MES_neftel2019 <- file.path(path_data, "env_data_MES_sctransform_neftel2019.RData")
path_NPC_neftel2019 <- file.path(path_data, "env_data_NPC_sctransform_neftel2019.RData")
path_OPC_neftel2019 <- file.path(path_data, "env_data_OPC_sctransform_neftel2019.RData")

# Set seed for example
set.seed(214)

load(path_general_env_data_neftel2019)
min_number_of_clusters <- 3
max_number_of_clusters <- 3
penalization_lambda <- c(0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5)
rm(d)  # If loaded remove it since it uses memory
gc()  # Force clean memory
target_gene_cluster_vector <- seq(min_number_of_clusters, max_number_of_clusters)

for (i_cell_cluster in c(3)) {
  if(i_cell_cluster==1){
    load(path_AC_neftel2019)
    current_cell_cluster <- cell_cluster_AC
    rm(cell_cluster_AC)
  }else if(i_cell_cluster==2){
    load(path_MES_neftel2019)
    current_cell_cluster <- cell_cluster_MES
    rm(cell_cluster_MES)
  }else if(i_cell_cluster==3){
    load(path_NPC_neftel2019)
    current_cell_cluster <- cell_cluster_NPC
    rm(cell_cluster_NPC)
  }else if(i_cell_cluster==4){
    load(path_OPC_neftel2019)
    current_cell_cluster <- cell_cluster_OPC
    rm(cell_cluster_OPC)
  }

  # Run screg with a bunch of different cluster setings
  results <- vector(mode = "list", length = length(target_gene_cluster_vector))

  for (i_n_target_genes_clusters in seq(length(target_gene_cluster_vector))) {
    n_target_genes_clusters <- target_gene_cluster_vector[i_n_target_genes_clusters]
    print(paste("Cell cluster", i_cell_cluster, "Number of target gene clusters",n_target_genes_clusters), quote=FALSE)

    results[[i_n_target_genes_clusters]] <- scregclust::scregclust(
      expression = current_cell_cluster,  # p rows of genes, n columns of cells
      genesymbols = rownames(current_cell_cluster),  # Gene row numbers
      is_regulator = is_regulator, # inverse_which(indices = ind_reggenes, output_length = n_regulator_genes + n_target_genes),  # Vector indicating which genes are regulators
      n_cl = n_target_genes_clusters,
      penalization = penalization_lambda,
      noise_threshold = 0.000001,
      verbose = TRUE,
      n_cycles = 18,
      compute_silhouette = TRUE,
      center=FALSE
    )
  }

  saveRDS(results, file.path(path_data, paste0("scregResultsNeftel2019_lambda_", paste0(penalization_lambda, collapse="-"),
                                               "_cellCluster_", i_cell_cluster,
                                               "_nRegulatorGenes_", n_regulator_genes,
                                               "_nTargetGenes_", n_target_genes,
                                               "_nCells_", ncol(current_cell_cluster),
                                               "_minNCluster_", min_number_of_clusters,
                                                "_maxNCluster", max_number_of_clusters,".rds")))
  rm(results, current_cell_cluster)
  gc()  # Force clean memory
}

