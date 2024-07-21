#!/usr/bin/Rscript
# This script is for producing plots and informaiton for multiple runs of scregclust
# Multiple runs are: n_cl varies and lambda varies
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

# ------------------ LOAD
# scregclust::plot_silhouettes(list_of_fits = results,
#                              penalization = penalization_lambda)
#
# scregclust::plot_cluster_count_helper(list_of_fits = results, penalization = penalization_lambda)
min_number_of_clusters <- 2
max_number_of_clusters <- 15
target_gene_cluster_vector <- seq(min_number_of_clusters, max_number_of_clusters)
penalization_lambda <- c(0.1, 0.2, 0.3, 0.5)
path_data <- "C:\\Users\\Sebastian\\repos\\biclust\\demo\\neftel_experiment"
results1 <- readRDS(file.path(path_data, "scregResultsNeftel2019_lambda_0.1-0.2-0.3-0.5_cellCluster_1_nRegulatorGenes_208_nTargetGenes_1640_nCells_4872_minNCluster_2_maxNCluster15.rds"))
results2 <- readRDS(file.path(path_data, "scregResultsNeftel2019_lambda_0.1-0.2-0.3-0.5_cellCluster_2_nRegulatorGenes_208_nTargetGenes_1640_nCells_2923_minNCluster_2_maxNCluster15.rds"))
results3 <- readRDS(file.path(path_data, "scregResultsNeftel2019_lambda_0.1-0.2-0.3-0.5_cellCluster_3_nRegulatorGenes_208_nTargetGenes_1640_nCells_2810_minNCluster_2_maxNCluster15.rds"))
results4 <- readRDS(file.path(path_data, "scregResultsNeftel2019_lambda_0.1-0.2-0.3-0.5_cellCluster_4_nRegulatorGenes_208_nTargetGenes_1640_nCells_2093_minNCluster_2_maxNCluster15.rds"))
names <- list("AC", "MES", "NPC", "OPC")
all_results <- list(results1, results2, results3, results4)



# ------------------ PLOTS

plot_things <- function() {
    tryCatch(
        {
             p <- scregclust::plot_cluster_count_helper(list_of_fits = all_results[[i]], penalization = cp) + ggtitle(paste(names[i], "lambda", cp)) + coord_cartesian(ylim = c(0, 1))
            # p <- scregclust::plot_silhouettes(list_of_fits = results[[1]], penalization = penalization_lambda)
            print(p)
        },
        error = function(cond) {
            NULL
        },
        warning = function(cond) {
            NULL
        },
        finally = {
            print("done")
        }
    )
}

plot_things2 <- function() {
    tryCatch(
        {
          p <- scregclust:::plot.scregclust(all_results[[i]][[ii]]) + ggtitle(paste(names[i], "n clusters", target_gene_cluster_vector[ii])) + coord_cartesian(ylim = c(0, 1))
          print(p)
        },
        error = function(cond) {
            NULL
        },
        warning = function(cond) {
            NULL
        },
        finally = {
            print("done")
        }
    )
}

for (i in seq(4)) {
  for (cp in penalization_lambda) {
    # print(str(as.numeric(cp)))
    print(paste(names[i], "lambda", cp), quote = FALSE)
    plot_things()

    # p <- scregclust::plot_cluster_count_helper(list_of_fits = all_results[[i]], penalization = cp) + ggtitle(paste(names[i], "lambda", cp))
    # # p <- scregclust::plot_silhouettes(list_of_fits = results[[1]], penalization = penalization_lambda)
    # print(p)
  }
}

for (i in seq(4)) {
  for (ii in seq(14)) {
    plot_things2()
  }
}






# -------------- GET INFO
# stats
# str(results1[[1]]$results[[1]]$output[[1]])
# colSums((results1[[1]]$results[[1]]$output[[1]]$reg_table)!=0)
# colSums(!is.na(results1[[1]]$results[[1]]$output[[1]]$importance))
# colSums((results1[[1]]$results[[1]]$output[[1]]$models))
#
# beta_logical <- results1[[1]]$results[[1]]$output[[1]]$models
# coeffs_old_style <- vector(mode = "list", length = ncol(beta_logical))
# coeffs_new_style <- results1[[1]]$results[[1]]$output[[1]]$coeffs
# for(i_target_gene_cluster in seq(length(coeffs_old_style))){
#   temp_coeffs <- matrix(0, nrow=nrow(beta_logical), ncol=ncol(coeffs_new_style[[i_target_gene_cluster]]))
#   current_beta_logical <- beta_logical[,i_target_gene_cluster]  # Reduces to a vector in R
#   temp_coeffs[current_beta_logical,] <- coeffs_new_style[[i_target_gene_cluster]]
#   coeffs_old_style[[i_target_gene_cluster]] <- temp_coeffs
# }
# rm(beta_logical, coeffs_new_style, temp_coeffs, current_beta_logical)
#



# Gets all the interesting info from a list of scregclust results
get_info <- function(res) {
  n_cluster_parameters <- length(res)
  n_sum_of_all_lambdas <- 0
  for (i_cluster in seq(n_cluster_parameters)) {
    n_sum_of_all_lambdas <- n_sum_of_all_lambdas + length(res[[i_cluster]]$results)
  }
  n_combos <- n_sum_of_all_lambdas

  n_cl <- vector(length = n_combos)
  penalization_lambdas <- vector(length = n_combos)
  r2 <- vector(length = n_combos)
  silhouette <- vector(length = n_combos)
  converged <- vector(length = n_combos)
  n_lost_clusters <- vector(length = n_combos)
  mean_regulators_in_non_empty_clusters <- vector(length = n_combos)
  n_noise_cluster <- vector(length = n_combos)

  regulator_index <- c()
  non_empty_target_gene_clusters <- 0
  i <- 0
  for (i_cluster in seq(n_cluster_parameters)) {
    n_lambdas <- length(res[[i_cluster]]$results)
    for (i_penalization in seq(n_lambdas)) {
      i <- i + 1
      n_cl[i] <- res[[i_cluster]]$results[[i_penalization]]$n_cl
      temp_r2 <- res[[i_cluster]]$results[[i_penalization]]$output[[1]]$r2_cluster
      r2[i] <- mean(temp_r2, na.rm = TRUE)
      temp_silhouette <- res[[i_cluster]]$results[[i_penalization]]$output[[1]]$silhouette
      silhouette[i] <- mean(temp_silhouette[is.finite(temp_silhouette)], na.rm = TRUE,)

      # Regulators in target gene cluster 1
      n_reg_genes <- vector(length = n_cl[i])
      for (i_target_gene_cluster in seq(n_cl[i])) {
        creg <- (res[[i_cluster]]$results[[i_penalization]]$output[[1]]$reg_table[[i_target_gene_cluster]])
        n_reg_genes[i_target_gene_cluster] <- sum(creg != 0, na.rm = TRUE)
        if (!all(is.na(creg))) {
          non_empty_target_gene_clusters <- non_empty_target_gene_clusters + 1
          regulator_index <- c(regulator_index, which(creg != 0))
        }
      }
      n_reg_genes_string <- paste(n_reg_genes, collapse = '-')
      n_actual_clusters <- sum(n_reg_genes != 0)
      mean_regulators_in_non_empty_clusters[i] <- mean(n_reg_genes[n_reg_genes != 0])

      # Target genes in gene cluster 1
      # sum(res[[1]]$results[[4]]$output[[1]]$cluster==1,na.rm=TRUE)
      n_noise_cluster[i] <- sum(res[[i_cluster]]$results[[i_penalization]]$output[[1]]$cluster == -1, na.rm = TRUE)
      n_lost_clusters[i] <- n_cl[i] - n_actual_clusters

      penalization_lambdas[i] <- res[[i_cluster]][[1]][[i_penalization]]
      converged[i] <- res[[i_cluster]]$results[[i_penalization]]$converged

      print(paste(n_cl[i],
                  penalization_lambdas[i],
                  r2[i],
                  silhouette[i],
                  converged[i],
                  "\"", n_reg_genes_string, "\"",
                  n_actual_clusters,
                  mean_regulators_in_non_empty_clusters[i],
                  n_noise_cluster[i],
                  n_lost_clusters[i]),
            quote = FALSE)
    }
  }

  df <- data.frame(n_cl,
                   penalization_lambdas,
                   r2,
                   silhouette,
                   converged,
                   n_lost_clusters,
                   mean_regulators_in_non_empty_clusters,
                   n_noise_cluster)
  df <- subset(df, converged == TRUE)
  df['rank'] <- rank((rank(-df['r2']) +
    rank(-df['silhouette']) +
    rank(-df['mean_regulators_in_non_empty_clusters']) +
    rank(df['n_lost_clusters']) +
    rank(df['n_noise_cluster'])))

  df <- df[order(df$rank),]
  rownames(df) <- NULL
  return(list("df" = df,
              "regulator_index" = regulator_index,
              "non_empty_target_gene_clusters" = non_empty_target_gene_clusters))
}

names <- list("AC", "MES", "NPC", "OPC")


current_results <- results3

RES <- get_info(res = current_results)
df <- RES['df']
regulator_index <- RES['regulator_index']$regulator_index
non_empty_target_gene_clusters <- RES$non_empty_target_gene_clusters
write.table(df, file = "", sep = ";", row.names = FALSE, col.names = TRUE, quote = FALSE)
regnames <- current_results[[1]]$results[[1]]$genesymbols
n_regulator_genes <- length(current_results[[1]]$results[[1]]$output[[1]]$reg_table[[1]])
names_of_regulator_genes <- regnames[(length(regnames)-n_regulator_genes+1):length(regnames)]
barplot(table(regulator_index)/non_empty_target_gene_clusters, ylim=c(0, 1))
d <- table(regulator_index)/non_empty_target_gene_clusters
d <- d[d>0.2]
d <- as.data.frame(d)
d['regulator_name'] <- names_of_regulator_genes[which(d>0.2)]


#
# for (i_cell_cluster in seq(n_cell_clusters)) {
#   if(i_cell_cluster==1){
#     load(path_AC_neftel2019)
#     current_cell_cluster <- cell_cluster_AC
#     print(dim(current_cell_cluster))
#     rm(cell_cluster_AC)
#   }else if(i_cell_cluster==2){
#     load(path_MES_neftel2019)
#     current_cell_cluster <- cell_cluster_MES
#      print(dim(current_cell_cluster))
#     rm(cell_cluster_MES)
#   }else if(i_cell_cluster==3){
#     load(path_NPC_neftel2019)
#     current_cell_cluster <- cell_cluster_NPC
#      print(dim(current_cell_cluster))
#     rm(cell_cluster_NPC)
#   }else if(i_cell_cluster==4){
#     load(path_OPC_neftel2019)
#     current_cell_cluster <- cell_cluster_OPC
#      print(dim(current_cell_cluster))
#     rm(cell_cluster_OPC)
#   }
# }
