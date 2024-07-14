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

# Ensemble gene IDs (like "ENSG00000269696") to standard gene names (like "A1BG"),
ensembl_to_normal <- function(data_matrix, ensembl_gene_names) {
  geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys = ensembl_gene_names, keytype = "GENEID", columns = c("SYMBOL", "GENEID"))
  select_these <- geneIDs1[, 2] # Since we probably lost genes in the translation
  normal_gene_names <- geneIDs1[, 1]

  if (nrow(data_matrix) == ncol(data_matrix)) {
    stop("Can't determine which dimension is genes/cells.")
  }

  if (nrow(data_matrix) == length(ensembl_gene_names)) {
    data_matrix <- data_matrix[select_these,]
  }else if (ncol(data_matrix) == length(ensembl_gene_names)) {
    data_matrix <- data_matrix[, select_these]
  }else(
    stop("The input ensembl vector doesn't match any dimension of data_matrix")
  )

  return(list("data_matrix" = data_matrix, "normal_gene_names" = normal_gene_names))
}

give_me_regulators <- function(normal_gene_names) {
  fake_matrix <- matrix(0, nrow = length(normal_gene_names), ncol = 1)  # This is fast
  rownames(fake_matrix) <- normal_gene_names  # as.vector(d@assays$SCT@var.features)
  out <- scregclust::scregclust_format(fake_matrix)  # Needs to be a matrix to use this
  is_regulator <- out[['is_regulator']]
  rm(fake_matrix, out)
  return(is_regulator)
}

if (file.exists(path_env_data_neftel2019) && file.exists(path_general_env_data_neftel2019)) {
  load(path_env_data_neftel2019)
}else {
  # Read data
  d <- readRDS(path_coreGBmap)

  # Keep the correct type of cells
  logical_keep_neftel_cells <- as.character(d@meta.data$author)=="Neftel2019"
  cell_types <- d@meta.data$annotation_level_3
  keep_cells <- cell_types == 'AC-like' | cell_types == 'MES-like' | cell_types == 'NPC-like' | cell_types == 'OPC-like'
  keep_cells <- logical_keep_neftel_cells & keep_cells
  rm(logical_keep_neftel_cells)

  cell_types <- factor(cell_types[keep_cells])
  d <- d[, keep_cells]
  rm(keep_cells)

  # Remove genes which are 'constant', aka no variance, just one number. Scregclust doesn't work on those.
  non_constant_ind <- which(apply(d@assays$RNA$data, 1, sd) != 0)
  d <- d[non_constant_ind,]
  rm(non_constant_ind)

  # Reduce the number of genes with sctransform. It converts gene names to ensembl gene names.
  d <- Seurat::SCTransform(d, variable.features.n = 2000)

  d <- d@assays$SCT$scale.data  # rows are genes, cols are cells

  # Remove all genes with less non zero numbers than 200
  d <- d[, apply(d, MARGIN = 2, function(x) sum(x != 0)) > 200]

  # Remove all genes/rows that don't correlate with other rows more than 0.1
  cor_matrix <- abs(cor(t(d)))  # For correlation we want cols as genes
  diag(cor_matrix) <- 0
  threshold <- 0.1
  keep_rows <- apply(cor_matrix, 1, max) > threshold
  d <- d[keep_rows,]
  rm(keep_rows, cor_matrix)

  # Ensembl gene IDs (like "ENSG00000269696") to standard gene names (like "A1BG"),
  geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= rownames(d), keytype = "GENEID", columns = c("SYMBOL","GENEID"))
  select_these <- geneIDs1[,2] # Since we probably lost genes in the translation
  new_names <- geneIDs1[,1]
  d <- d[select_these,]
  rownames(d) <- new_names
  rm(select_these)

  # Find out which genes/rows are regulators.
  is_regulator <- give_me_regulators(normal_gene_names = rownames(d))

  # Put gene/rows in order. This puts the regulator genes at the end.
  d <- d[order(is_regulator),]
  is_regulator <- is_regulator[order(is_regulator)]

  # Put cells/columns in order.
  d <- d[, order(cell_types)]
  cell_types <- cell_types[order(cell_types)]

  n_regulator_genes <- sum(is_regulator)
  n_target_genes <- nrow(d) - n_regulator_genes


  # Now when we have put the cells in order we just need to count the cells
  # in each cell cluster. Then it's easy to make a vector with the true_cluster_allocation
  # further down in code
  n_cells_cell_cluster_1 <- sum(cell_types == 'AC-like')
  n_cells_cell_cluster_2 <- sum(cell_types == 'MES-like')
  n_cells_cell_cluster_3 <- sum(cell_types == 'NPC-like')
  n_cells_cell_cluster_4 <- sum(cell_types == 'OPC-like')

  # Print some stats
  print(object.size(d), units = 'MB', standard = 'SI')
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
  true_cluster_allocation <- rep(1:n_cell_clusters, times = n_cells_in_each_cluster)
  rm(n_cells_in_each_cluster, n_cells_cell_cluster_1, n_cells_cell_cluster_2, n_cells_cell_cluster_3, n_cells_cell_cluster_4)

  ind_targetgenes <- which(c(rep(1, n_target_genes), rep(0, n_regulator_genes)) == 1)
  ind_reggenes <- which(c(rep(0, n_target_genes), rep(1, n_regulator_genes)) == 1)

  gc()  # Force clean memory


  save(d,
       is_regulator,
       n_cell_clusters,
       n_regulator_genes,
       n_target_genes,
       ind_reggenes,
       ind_targetgenes,
       true_cluster_allocation,
       file = path_env_data_neftel2019)

   save(is_regulator,
       n_cell_clusters,
       n_regulator_genes,
       n_target_genes,
       ind_reggenes,
       ind_targetgenes,
       true_cluster_allocation,
       file = path_general_env_data_neftel2019)
}

# Create seperate datasets for each cell state
if (!file.exists(path_AC_neftel2019) ||
  !file.exists(path_MES_neftel2019) ||
  !file.exists(path_NPC_neftel2019) ||
  !file.exists(path_OPC_neftel2019)) {
  print(object.size(d), units = 'MB', standard = 'SI')

  load(path_env_data_neftel2019)

  # Save one variable per cell cluster to save on memory use
  cell_cluster_AC <- d[, true_cluster_allocation == 1]
  print(dim(cell_cluster_AC))
  print(object.size(cell_cluster_AC), units = 'MB', standard = 'SI')
  save(cell_cluster_AC, file = path_AC_neftel2019)
  rm(cell_cluster_AC)
  gc()  # Force clean memory

  cell_cluster_MES <- d[, true_cluster_allocation == 2]
  print(dim(cell_cluster_MES))
  print(object.size(cell_cluster_MES), units = 'MB', standard = 'SI')
  save(cell_cluster_MES, file = path_MES_neftel2019)
  rm(cell_cluster_MES)
  gc()  # Force clean memory

  cell_cluster_NPC <- d[, true_cluster_allocation == 3]
  print(dim(cell_cluster_NPC))
  print(object.size(cell_cluster_NPC), units = 'MB', standard = 'SI')
  save(cell_cluster_NPC, file = path_NPC_neftel2019)
  rm(cell_cluster_NPC)
  gc()  # Force clean memory

  cell_cluster_OPC <- d[, true_cluster_allocation == 4]
  print(dim(cell_cluster_OPC))
  print(object.size(cell_cluster_OPC), units = 'MB', standard = 'SI')
  save(cell_cluster_OPC, file = path_OPC_neftel2019)
  rm(cell_cluster_OPC)
  gc()  # Force clean memory
}

# ----------------------------
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


# TODO: Not updated yet
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
