#!/usr/bin/Rscript
# This script is for producing plots and informaiton for multiple runs of scregclust
# Multiple runs are: n_modules varies and lambda varies
Sys.setlocale("LC_CTYPE", "en_US.UTF-8")

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

# Set seed for example
set.seed(214)

# Load data
files <- list.files(
  path = local_data, # Current directory, change if needed
  pattern = "^scregResults_pbmc10k_lambda_0",
  full.names = TRUE
)

names <- paste0("Cluster_", seq(length(files)))

# ------------------ PLOTS

plot_things <- function() {
    tryCatch(
        {
             p <- scregclust::plot_module_count_helper(list_of_fits = all_results[[i]], penalization = cp) + ggtitle(paste(names[i], "lambda", cp)) + coord_cartesian(ylim = c(0, 1))
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
          p <- scregclust:::plot.scregclust(all_results[[i]][[ii]]) + ggtitle(paste(names[i], "n modules", target_gene_module_vector[ii])) + coord_cartesian(ylim = c(0, 1))
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

    # p <- scregclust::plot_module_count_helper(list_of_fits = all_results[[i]], penalization = cp) + ggtitle(paste(names[i], "lambda", cp))
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
# for(i_target_gene_module in seq(length(coeffs_old_style))){
#   temp_coeffs <- matrix(0, nrow=nrow(beta_logical), ncol=ncol(coeffs_new_style[[i_target_gene_module]]))
#   current_beta_logical <- beta_logical[,i_target_gene_module]  # Reduces to a vector in R
#   temp_coeffs[current_beta_logical,] <- coeffs_new_style[[i_target_gene_module]]
#   coeffs_old_style[[i_target_gene_module]] <- temp_coeffs
# }
# rm(beta_logical, coeffs_new_style, temp_coeffs, current_beta_logical)
#



# Gets all the interesting info from a list of scregclust results
get_info <- function(res) {
  n_module_parameters <- length(res)
  n_sum_of_all_lambdas <- 0
  for (i_module in seq(n_module_parameters)) {
    n_sum_of_all_lambdas <- n_sum_of_all_lambdas + length(res[[i_module]])
  }
  n_combos <- n_sum_of_all_lambdas

  n_modules <- vector(length = n_combos)
  penalization_lambdas <- vector(length = n_combos)
  r2 <- vector(length = n_combos)
  silhouette <- vector(length = n_combos)
  converged <- vector(length = n_combos)
  n_lost_modules <- vector(length = n_combos)
  mean_regulators_in_non_empty_modules <- vector(length = n_combos)
  n_noise_modules <- vector(length = n_combos)

  regulator_index <- c()
  non_empty_target_gene_modules <- 0
  i <- 0
  for (i_module in seq(n_module_parameters)) {
    n_lambdas <- length(res[[i_module]])
    for (i_penalization in seq(n_lambdas)) {
      i <- i + 1
      if(is.na(current_results[[i_module]][[i_penalization]][[1]])){
        next
      }
      if(is.null(res[[i_module]][[i_penalization]]$results[[1]])){
        next
      }

      n_modules[i] <- res[[i_module]][[i_penalization]]$results[[1]]$n_modules
      temp_r2 <- res[[i_module]][[i_penalization]]$results[[1]]$output[[1]]$r2_module
      r2[i] <- mean(temp_r2, na.rm = TRUE)
      temp_silhouette <- res[[i_module]][[i_penalization]]$results[[1]]$output[[1]]$silhouette
      silhouette[i] <- mean(temp_silhouette[is.finite(temp_silhouette)], na.rm = TRUE,)

      # Regulators in target gene module 1
      n_reg_genes <- vector(length = n_modules[i])
      for (i_target_gene_module in seq(n_modules[i])) {
        creg <- (res[[i_module]][[i_penalization]]$results[[1]]$output[[1]]$reg_table[[i_target_gene_module]])
        n_reg_genes[i_target_gene_module] <- sum(creg != 0, na.rm = TRUE)
        if (!all(is.na(creg))) {
          non_empty_target_gene_modules <- non_empty_target_gene_modules + 1
          regulator_index <- c(regulator_index, which(creg != 0))
        }
      }
      n_reg_genes_string <- paste(n_reg_genes, collapse = '-')
      n_actual_modules <- sum(n_reg_genes != 0)
      mean_regulators_in_non_empty_modules[i] <- mean(n_reg_genes[n_reg_genes != 0])

      # Target genes in gene module 1
      # sum(res[[1]]$results[[4]]$output[[1]]$module==1,na.rm=TRUE)
      n_noise_modules[i] <- sum(res[[i_module]][[i_penalization]]$results[[1]]$output[[1]]$module == -1, na.rm = TRUE)
      n_lost_modules[i] <- n_modules[i] - n_actual_modules

      penalization_lambdas[i] <- res[[i_module]][[i_penalization]][[1]]
      converged[i] <- res[[i_module]][[i_penalization]]$results[[1]]$converged

      # print(paste(n_modules[i],
      #             penalization_lambdas[i],
      #             r2[i],
      #             silhouette[i],
      #             converged[i],
      #             "\"", n_reg_genes_string, "\"",
      #             n_actual_modules,
      #             mean_regulators_in_non_empty_modules[i],
      #             n_noise_modules[i],
      #             n_lost_modules[i]),
      #       quote = FALSE)
    }
  }

  df <- data.frame(n_modules,
                   penalization_lambdas,
                   r2,
                   silhouette,
                   converged,
                   n_lost_modules,
                   mean_regulators_in_non_empty_modules,
                   n_noise_modules)
  df <- subset(df, converged == TRUE)
  df['rank'] <- rank((rank(-df['r2']) +
    rank(-df['silhouette']) +
    rank(-df['mean_regulators_in_non_empty_modules']) +
    rank(df['n_lost_modules']) +
    rank(df['n_noise_modules'])))

  df <- df[order(df$rank),]
  rownames(df) <- NULL
  return(list("df" = df,
              "regulator_index" = regulator_index,
              "non_empty_target_gene_modules" = non_empty_target_gene_modules))
}


for(cf in files){
  current_results <- readRDS(cf)
  RES <- get_info(res = current_results)
  df <- RES['df']
  regulator_index <- RES['regulator_index']$regulator_index
  non_empty_target_gene_modules <- RES$non_empty_target_gene_modules
  write.table(df, file = "", sep = ";", row.names = FALSE, col.names = TRUE, quote = FALSE)
}



regnames <- NA
for (i_nModules in seq(length(current_results))) {

  for (i_nLambda in seq(length(current_results[[i_nModules]]))) {
    if (any(is.na(current_results[[i_nModules]][[i_nLambda]]))) {
      next
    }else {
      regnames <- current_results[[i_nModules]][[i_nLambda]]$results[[1]]$genesymbols
      n_regulator_genes <- length(current_results[[i_nModules]][[i_nLambda]]$results[[1]]$output[[1]]$reg_table[[1]])
      break
    }
  }
  if(length(regnames)>1){
    break
  }
}

names_of_regulator_genes <- regnames[(length(regnames)-n_regulator_genes+1):length(regnames)]
barplot(table(regulator_index)/non_empty_target_gene_modules, ylim=c(0, 1))
temp_data <- table(regulator_index)/non_empty_target_gene_modules
temp_data <- temp_data[temp_data>0]
temp_names <- names_of_regulator_genes[which(temp_data>0)]
temp_data <- as.data.frame(temp_data)
temp_data['regulator_name'] <- temp_names
write.table(temp_data, file = "", sep = ";", row.names = FALSE, col.names = TRUE, quote = FALSE)

#
# for (i_cell_module in seq(n_cell_modules)) {
#   if(i_cell_module==1){
#     load(path_AC_neftel2019)
#     current_cell_module <- cell_module_AC
#     print(dim(current_cell_module))
#     rm(cell_module_AC)
#   }else if(i_cell_module==2){
#     load(path_MES_neftel2019)
#     current_cell_module <- cell_module_MES
#      print(dim(current_cell_module))
#     rm(cell_module_MES)
#   }else if(i_cell_module==3){
#     load(path_NPC_neftel2019)
#     current_cell_module <- cell_module_NPC
#      print(dim(current_cell_module))
#     rm(cell_module_NPC)
#   }else if(i_cell_module==4){
#     load(path_OPC_neftel2019)
#     current_cell_module <- cell_module_OPC
#      print(dim(current_cell_module))
#     rm(cell_module_OPC)
#   }
# }
