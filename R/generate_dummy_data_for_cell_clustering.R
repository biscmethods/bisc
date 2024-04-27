R_path <- here::here("R")
source(paste0(R_path, "/generate_dummy_data_for_scregclust.R"))
library(plyr)

# # dev
# n_cell_clusters <- 2
# n_target_gene_clusters <- c(2, 3)  # Number of target gene clusters in each cell cluster
# n_target_genes <- 30
# n_regulator_genes <- 100
# n_cells <- c(10000, 10000)
# regulator_means <- c(5, 1)  # For generating dummy data, regulator mean in each cell cluster
# coefficient_means <- list(c(1, 1.2), c(5, 5.5, 6))  # For generating dummy data, coefficient means in each cell cluster
# coefficient_sds <- list(c(0.1, 0.1), c(0.1, 0.1, 0.1))
# true_cluster_allocation <- rep(1:n_cell_clusters, times = n_cells)
# disturbed_fraction = 0.1


generate_dummy_data_for_cell_clustering <- function(
  n_cell_clusters = 3,
  n_target_gene_clusters = c(3, 4, 5),  # Number of target gene clusters in each cell cluster
  n_target_genes = 40,
  n_regulator_genes = 20,
  n_cells = c(1000, 5000, 10000),
  regulator_means = c(1, 2, 3),  # For generating dummy data, regulator mean in each cell cluster
  coefficient_means = list(c(1, 20, 30), c(1, 2, 3, 4), c(1, 2, 3, 4, 5)),  # For generating dummy data, coefficient means in each cell cluster
  coefficient_sds = list(c(0.1, 0.1, 0.1), c(0.1, 0.1, 0.1, 0.1), c(0.1, 0.1, 0.1, 0.1, 0.1)),
  disturbed_fraction = 0.1,  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
  plot_stuff = FALSE
) {

  # Set variables ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  true_cluster_allocation <- rep(1:n_cell_clusters, times = n_cells)
  total_n_cells <- sum(n_cells)

  # Generate dummy data for each cell cluster that we want ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  dummy_data <- vector(mode = "list", length = n_cell_clusters)
  betas <- dummy_data
  for (i_cluster in 1:n_cell_clusters) {
    print(paste(" Generating data for cell cluster", i_cluster), quote = FALSE)
    dummy_data[[i_cluster]] <- generate_dummy_data_for_scregclust(n_target_genes,
                                                                  n_regulator_genes,
                                                                  n_cells = n_cells[i_cluster],
                                                                  n_target_gene_clusters = n_target_gene_clusters[i_cluster],
                                                                  regulator_mean = regulator_means[i_cluster],
                                                                  coefficient_mean = coefficient_means[[i_cluster]],
                                                                  coefficient_sd = coefficient_sds[[i_cluster]])
    # list of list of models. with dimension regulators x targets
    # Top level list is each cell cluster
    # in each of those elements is the model generating each target gene cluster in that cell cluster
    betas[[i_cluster]] <- dummy_data[[i_cluster]]$B
  }



  # Create Z_r and Z_t from dummy data
  Z_t <- dummy_data[[1]]$Z_t  # cells x genes
  Z_r <- dummy_data[[1]]$Z_r  # cells x genes
  if (n_cell_clusters > 1) {
    for (i_cluster in 2:n_cell_clusters) {
      Z_t <- rbind(Z_t, dummy_data[[i_cluster]]$Z_t)
      Z_r <- rbind(Z_r, dummy_data[[i_cluster]]$Z_r)
    }
  }

  #
  # regulator_genes_index <- (1:(n_target_genes + n_regulator_genes) > n_target_genes) + 0

  # scregclust::scregclust(
  #   expression = rbind(t(Z_t), t(Z_r)),    #scRegClust wants this form
  #   genesymbols = 1:(n_target_genes + n_regulator_genes),               #gene row numbers
  #   is_regulator = regulator_genes_index, #vector indicating which genes are regulators
  #   n_cl = n_target_gene_clusters[i_cluster],
  #   penalization = 0.001,
  #   verbose = FALSE
  # ) -> scRegOut

  # str(scRegOut)
  # str(scRegOut$results)
  # str(scRegOut$results[[1]])
  # str(scRegOut$results[[1]]$output)
  # scRegOut$results[[1]]$output[[1]]$coeffs
  # scRegOut$results[[1]]$output[[1]]$models
  # scRegOut$results[[1]]$output[[1]]$cluster[!regulator_genes_index]

  dat <- cbind(Z_t, Z_r)

  # Split into train and test data for cell clustering
  # Skip for now
  # cell_data_split    <- sample(c(1,2), nrow(Z_t), replace = T)
  # train_indices      <- which(cell_data_split == 1)
  # train_dat          <- dat[,train_indices]

  # Get initial cell clustering
  initial_cell_clust <- kmeans(dat, n_cell_clusters)$cluster

  # initial_cell_clust <- sample(1:n_cell_clusters, n_cells, replace = T)


  # Kod för att flytta 1% av cellerna i varje kluster till ett annat kluster.
  disturbed_initial_cell_clust <- true_cluster_allocation
  if (disturbed_fraction > 0) {
    for (i_cluster in 1:n_cell_clusters) {
      indexes_of_cluster <- which(true_cluster_allocation == i_cluster)
      some_of_those_indexes <- sample(indexes_of_cluster, size = as.integer(length(indexes_of_cluster) * disturbed_fraction), replace = F)
      disturbed_initial_cell_clust[some_of_those_indexes] <- sample((1:n_cell_clusters)[-i_cluster], size = length(some_of_those_indexes), replace = T)
    }
  }



  if(plot_stuff){
    #plot regulators and targets, color by true cell cluster

    demo_path <- here::here("demo")
    R_path <- here::here("R")
    output_path <- demo_path

    library(ggmulti)
    library(ggplot2)
    library(dplyr)
    library(ggfortify)


    # plot regulators and targets
    plot_pca <- function(indata     = vis_dat$regulators,
                         cluster_allocations = vis_dat$true_cell_cluster_allocation,
                         plot_title = "Principal components of regulator genes"
    ){

      library(ggmulti)
      library(ggplot2)
      library(dplyr)
      library(ggfortify)

      # Apply PCA
      pca_result <- prcomp(indata, scale. = TRUE)

      # Extract principal components
      pc_df <- as_tibble(pca_result$x[,1:2])

      # Add labels (if you have any)
      pc_df$label <- as.factor(cluster_allocations)

      # Plot the first two principal components
      plot_name <- "regulators_pca.png"

      # png(file.path(demo_path,plot_name))
      pc_df[sample(nrow(pc_df), 1000),] %>%
        ggplot( aes(x = PC1, y = PC2, colour = label)) +
        geom_point(alpha = 0.1) +
        theme_minimal() +
        labs(title = plot_title) +
        scale_color_discrete(name = "Cell cluster")
      # dev.off()
    }

    png(file.path(output_path, paste0("biclust_data_gen_regulators.png")),
        width = 1024, height = 480, units = "px")
    plot_pca(indata     = Z_r,
             cluster_allocations = true_cluster_allocation,
             plot_title = "Principal components of regulator genes") -> p
    print(p)
    dev.off()

    png(file.path(output_path, paste0("biclust_data_gen_targets.png")),
        width = 1024, height = 480, units = "px")
    plot_pca(indata     = Z_t,
             cluster_allocations = true_cluster_allocation,
             plot_title = "Principal components of target genes") -> p
    print(p)
    dev.off()


    # Visualise S to show cluster structures

    df_list <- lapply(seq_along(dummy_data), function(i) {
      df <- reshape::melt(dummy_data[[i]]$S)
      colnames(df) <- c("target_gene_cluster", "regulator", "value")
      df$plot <- paste("Cell cluster", i)

      return(df)
    })

    # Combine all data frames into one
    df_combined <- do.call(rbind, df_list)

    # Create a faceted plot
    ggplot(df_combined, aes(x = regulator, y = as.factor(target_gene_cluster), fill = value)) +
      geom_tile() +
      facet_wrap(~ plot, ncol = 2) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            plot.title = element_text(hjust = 0.5)) +
      scale_fill_gradient2(low = "#075AFF",
                           mid = "#FFFFCC",
                           high = "#FF0000") +
      xlab("Regulator gene") +
      ylab("Target gene cluster")  +
      ggtitle("Heatmap of S matrices") -> p


    png(file.path(output_path, paste0("biclust_data_gen_S.png")),
        width = 1024, height = 480, units = "px")
    print(p)
    dev.off()

    #visualise target gene cluster structure somehow?

    df_list <- lapply(seq_along(dummy_data), function(i) {
      df <- reshape::melt(dummy_data[[i]]$Pi)
      colnames(df) <- c("target_gene_cluster", "target_gene", "value")
      df$plot <- paste("Cell cluster", i)

      return(df)
    })

    # Combine all data frames into one
    df_combined <- do.call(rbind, df_list)

    ggplot(df_combined, aes(x = target_gene,
                            y = as.factor(target_gene_cluster),
                            fill = value)) +
      geom_tile(color = "black") +
      facet_wrap(~ plot, ncol = 2) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            plot.title = element_text(hjust = 0.5)) +
      scale_fill_gradient2(low = "#075AFF",
                           mid = "#FFFFCC",
                           high = "#FF0000") +
      xlab("Target gene") +
      ylab("Target gene cluster")  +
      ggtitle("Heatmap of Pi matrices") -> p

    png(file.path(output_path, paste0("biclust_data_gen_Pp.png")),
          width = 1024, height = 480, units = "px")
    print(p)
    dev.off()


    #plot betas (as heatmaps of matrix??? idk)


    df_list <- lapply(seq_along(dummy_data), function(i) {
      df <- reshape::melt(dummy_data[[i]]$B)
      colnames(df) <- c("regulator_gene", "target_gene", "value")
      df$plot <- paste("Cell cluster", i)

      return(df)
    })

    # Combine all data frames into one
    df_combined <- do.call(rbind, df_list)

    ggplot(df_combined, aes(x = regulator_gene,
                            y = target_gene, fill = value)) +
      geom_tile() +
      facet_wrap(~ plot, ncol = 2) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            plot.title = element_text(hjust = 0.5)) -> p

    png(file.path(output_path, paste0("biclust_data_gen_B.png")),
        width = 1024, height = 480, units = "px")
    print(p)
    dev.off()

  }

  return(list(disturbed_initial_cell_clust = disturbed_initial_cell_clust,
              initial_cell_clust = initial_cell_clust,
              true_cell_clust = true_cluster_allocation,
              true_betas = betas,
              dat = dat))
}
