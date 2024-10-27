#!/usr/bin/Rscript
rm(list = ls())

library(here)  # To work with paths
library(callr)
# library(patchwork)
library(scregclust)
# To be able to run library("EnsDb.Hsapiens.v79") you need to at least:
# BiocManager::install("GenomeInfoDb")
# BiocManager::install("SparseArray")
# BiocManager::install("EnsDb.Hsapiens.v79")
library("EnsDb.Hsapiens.v79")
sink()  # Because some of our scripts redirects output from scregclust to force it to be quite. This restores output.
gc()  # Force clean memory
Sys.setlocale("LC_CTYPE", "en_US.UTF-8")

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

# Function to monitor the output and stop if "apple" is detected
monitor_function <- function(function_to_run, ...) {
  cat("START\n")
  output <- character()

  # Start the process
  r_process <- callr::r_bg(func=function(function_to_run, ...) function_to_run(...),
                           stdout = "|",
                           stderr = "|",
                           supervise=TRUE,
                           args=list(function_to_run, ...))

  # Monitor the process output
  while (r_process$is_alive()) {
    # Read the standard output
    new_output <- r_process$read_output_lines()
    new_error <- r_process$read_error_lines()
    #cat(new_output)
    if (length(new_output) > 0) {
      new_output_merged <- stringr::str_trim(paste(new_output, collapse = "\n"))
      new_output_merged_clean <- gsub("[\r\n]", "", new_output_merged)
      if(!new_output_merged_clean==""){
        # new_output_merged <- stringi::stri_enc_toutf8(new_output_merged)
        # print(new_output_merged, quote=FALSE)  # Print the output to the console (optional)
        cat(new_output, sep = "\n")
        output <- c(output, new_output)  # Append new output to the variable


      }
    }

    #cat(new_output)
    if (length(new_error) > 0) {
      new_error_merged <- paste(new_error, collapse = "\n")
      new_error_merged_clean <- gsub("[\r\n]", "", new_error_merged)
      if (!new_error_merged_clean == "") {
        cat("ERROR: ")
        cat(new_error, sep = "\n")  # Print the output to the console (optional)
        output <- c(output, new_error)
      }
    }

    # Check if "apple" is in the new output
    if (any(grepl("Max iterations", new_output)) || any(grepl("Max iterations", new_error))) {
      r_process$kill()
      cat("Detected 'Max iterations'. Terminating process.\n")
      break
    }

    Sys.sleep(0.01)  # Small delay to prevent busy-waiting
  }

  cat("END\n")

  output_merged <- paste(output, collapse = "\n")
  cat(output_merged)

  # Capture any remaining output if process completed without printing "apple"
  if (!r_process$is_alive()) {
    result <- r_process$get_result()
  }

  # Convert output to a single string (if needed)
  # output <- paste(output, collapse = "\n")
  return(list("result"=result, "output"=output_merged))
}

# fruits <- function() {
#   cat("orange\n")
#   cat("banana\n")
#   cat("apple\n")
#   sleep(1)
#   cat("grape\n")
#   a <- "KNAPP"
#   return(a)
# }
#
# # Run the monitor function
# captured_output <- monitor_function(fruits)
# cat("Captured output:\n", captured_output)

#
# exit()


load(path_general_env_data_neftel2019)
min_number_of_clusters <- 15
max_number_of_clusters <- 20
penalization_lambda <- c(0.1,0.2,0.3,0.5)
rm(d)  # If loaded remove it since it uses memory
gc()  # Force clean memory
target_gene_cluster_vector <- seq(min_number_of_clusters, max_number_of_clusters)

for (i_cell_cluster in c(1,4,5)) {
  if (i_cell_cluster == 1) {
    load(path_AC_neftel2019)
    current_cell_cluster <- cell_cluster_AC
    rm(cell_cluster_AC)
  }else if (i_cell_cluster == 2) {
    load(path_MES_neftel2019)
    current_cell_cluster <- cell_cluster_MES
    rm(cell_cluster_MES)
  }else if (i_cell_cluster == 3) {
    load(path_NPC_neftel2019)
    current_cell_cluster <- cell_cluster_NPC
    rm(cell_cluster_NPC)
  }else if (i_cell_cluster == 4) {
    load(path_OPC_neftel2019)
    current_cell_cluster <- cell_cluster_OPC
    rm(cell_cluster_OPC)
  }else if (i_cell_cluster ==5){
    load(path_AC_neftel2019)
    load(path_OPC_neftel2019)
    current_cell_cluster <- cbind(cell_cluster_AC, cell_cluster_OPC)
  }

  # Run screg with a bunch of different cluster setings
  results <- vector(mode = "list", length = length(target_gene_cluster_vector))
  outputs <- vector(mode = "list", length = length(target_gene_cluster_vector))

  for (i_n_target_genes_clusters in seq(length(target_gene_cluster_vector))) {
    n_target_genes_clusters <- target_gene_cluster_vector[i_n_target_genes_clusters]
    print(paste("Cell cluster", i_cell_cluster, "Number of target gene clusters", n_target_genes_clusters), quote = FALSE)
    results[[i_n_target_genes_clusters]] <- vector(mode = "list", length = length(penalization_lambda))
    outputs[[i_n_target_genes_clusters]] <- vector(mode = "list", length = length(penalization_lambda))
    for (i_penalization_lambda in seq(length(penalization_lambda))) {
      current_penalization_lambda <- penalization_lambda[i_penalization_lambda]
      temp_output <- monitor_function(function_to_run = scregclust::scregclust,
                                      expression = current_cell_cluster,  # p rows of genes, n columns of cells
                                      genesymbols = rownames(current_cell_cluster),  # Gene row numbers
                                      is_regulator = is_regulator, # inverse_which(indices = ind_reggenes, output_length = n_regulator_genes + n_target_genes),  # Vector indicating which genes are regulators
                                      n_cl = n_target_genes_clusters,
                                      penalization = current_penalization_lambda,
                                      verbose = TRUE,
                                      n_cycles = 300,
                                      compute_silhouette = TRUE,
                                      center = FALSE)
      results[[i_n_target_genes_clusters]][[i_penalization_lambda]] <- temp_output[['result']]
      outputs[[i_n_target_genes_clusters]][[i_penalization_lambda]] <- temp_output[['output']]

    }
  }

  saveRDS(results, file.path(path_data, paste0("scregResultsNeftel2019_lambda_", paste0(penalization_lambda, collapse = "-"),
                                                           "_cellCluster_", i_cell_cluster,
                                                           "_nRegulatorGenes_", n_regulator_genes,
                                                           "_nTargetGenes_", n_target_genes,
                                                           "_nCells_", ncol(current_cell_cluster),
                                                           "_minNCluster_", min_number_of_clusters,
                                                           "_maxNCluster", max_number_of_clusters, ".rds")))
  saveRDS(outputs, file.path(path_data, paste0("scregResultsNeftel2019_lambda_printOutput_", paste0(penalization_lambda, collapse = "-"),
                                                           "_cellCluster_", i_cell_cluster,
                                                           "_nRegulatorGenes_", n_regulator_genes,
                                                           "_nTargetGenes_", n_target_genes,
                                                           "_nCells_", ncol(current_cell_cluster),
                                                           "_minNCluster_", min_number_of_clusters,
                                                           "_maxNCluster", max_number_of_clusters, ".rds")))

  rm(results, current_cell_cluster)
  gc()  # Force clean memory
}

