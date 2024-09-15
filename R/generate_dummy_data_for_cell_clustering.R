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
  n_cell_clusters = 2,
  n_target_gene_clusters = c(10, 7),  # Number of target gene clusters in each cell cluster
  n_target_genes = 2193,          #from vignette
  n_regulator_genes = 493,        # from vignette
  n_cells = c(6000, 6000),
  regulator_means = c(0, 0),  # For generating dummy data, regulator mean in each cell cluster
  coefficient_means = list(c(0.0417,
                             0.0343,
                             0.0576,
                             0.043 ,
                             0.0576,
                             0.0413,
                             0.0473,
                             0.0444,
                             0.0481,
                            -0.0139),
                           c(0.0404,
                             0.0519,
                             0.0545,
                             0.0915,
                             0.0663,
                             0.0512,
                            -0.0064
                            )
                           ),  # For generating dummy data, coefficient means in each cell cluster
  coefficient_sds = list(c(0.0556,
                           0.037,
                           0.0638,
                           0.0466,
                           0.0761,
                           0.0471,
                           0.0468,
                           0.0611,
                           0.0623,
                           0.0394
                           ),
                         c(0.0462,
                           0.0496,
                           0.0807,
                           0.1086,
                           0.071,
                           0.0716,
                           0.057
                           )
                         ),
  disturbed_fraction = 0.1,  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
  plot_stuff = TRUE,
  plot_suffix = "vignette",
  testing_penalization = c(0.1, 0.3) #vector of length n_cell_clusters
) {

  # Set variables ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  true_cluster_allocation <- rep(1:n_cell_clusters, times = n_cells)
  total_n_cells <- sum(n_cells)

  # Generate dummy data for each cell cluster that we want ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  dummy_data <- vector(mode = "list", length = n_cell_clusters)
  betas      <- dummy_data
  for (i_cluster in 1:n_cell_clusters) {
    print(paste(" Generating data for cell cluster", i_cluster), quote = FALSE)
    dummy_data[[i_cluster]] <- generate_dummy_data_for_scregclust(n_target_genes,
                                                                  n_regulator_genes,
                                                                  n_cells = n_cells[i_cluster],
                                                                  n_target_gene_clusters = n_target_gene_clusters[i_cluster],
                                                                  regulator_mean = regulator_means[i_cluster],
                                                                  coefficient_mean = coefficient_means[[i_cluster]],
                                                                  coefficient_sd = coefficient_sds[[i_cluster]],
                                                                  make_regulator_network_plot = F,
                                                                  plot_suffix = paste0(plot_suffix, "_", i_cluster),
                                                                  testing_penalization = testing_penalization[i_cluster]
                                                                  )
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

  # print("now jump into debug")
  # browser()
  ###########################################################################
  ##################### Create corresponding count data #####################
  ###########################################################################
  num_genes <- n_target_genes + n_regulator_genes
  num_cells <- sum(n_cells)

  theta <- runif(n = num_genes,
                 min = 0,
                 max = 1
  )  #dispersion per gene

  avg_counts_per_cell <-  30 / num_genes  # sensitivity of the sequencing machine isch

  counts <- matrix(data = 0, nrow = num_cells, ncol = num_genes)

  temp_dat <-dat/max(dat) #hack to make sensible counts for now

  for(cell in 1:num_cells){
    for(gene in 1:num_genes){
      counts[cell, gene] <- rnegbin(1,
                                    mu = avg_counts_per_cell * exp(temp_dat[cell, gene]),
                                    theta = theta[gene])
    }
  }

  ##########################################################################
  ###################### Do some other clustering stuff? ###################
  ##############################Consider removing ##########################

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

    cat("plotting stuff for the simulated cell cluestering data")

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
                         plot_title = "PCA plot"
    ){

      library(ggmulti)
      library(ggplot2)
      library(dplyr)
      library(ggfortify)
      library(Rtsne)

      # Apply PCA
      pca_result <- stats::prcomp(indata, scale. = TRUE)

      # Extract principal components
      pc_df <- as_tibble(pca_result$x[,1:2])

      # Add labels (if you have any)
      pc_df$label <- as.factor(cluster_allocations)

      pc_df[sample(nrow(pc_df), 1000),] %>%
        ggplot( aes(x = PC1, y = PC2, colour = label)) +
        geom_point(alpha = 0.4) +
        theme_minimal() +
        labs(title = plot_title) +
        scale_color_discrete(name = "Cell cluster")
      # dev.off()
    }

    png(file.path(output_path, paste0("biclust_data_gen_regulators_",plot_suffix,".png")),
        width = 1024, height = 480, units = "px")
    plot_pca(indata     = Z_r,
             cluster_allocations = true_cluster_allocation,
             plot_title = "Principal components of regulator genes") -> p
    print(p)
    dev.off()

    png(file.path(output_path, paste0("biclust_data_gen_targets_",plot_suffix,".png")),
        width = 1024, height = 480, units = "px")
    plot_pca(indata     = Z_t,
             cluster_allocations = true_cluster_allocation,
             plot_title = "Principal components of target genes") -> p
    print(p)
    dev.off()

    png(file.path(output_path, paste0("biclust_data_gen_complete_",plot_suffix,".png")),
        width = 1024, height = 480, units = "px")
    plot_pca(indata     = cbind(Z_t,Z_r),
             cluster_allocations = true_cluster_allocation,
             plot_title = "Principal components plot of simulated data") -> p
    print(p)
    dev.off()

    tsne_plot <- function(
      indata = dat,
      cluster_allocation = true_cluster_allocation,
      title = "tSNE plot"
    ){
      tsne_result <- Rtsne::Rtsne(
        indata, dims = 2,
        perplexity = 30,
        max_iter = 500,
        partial_pca = T,
        verbose = F
      )
      tsne_data   <- data.frame(X = tsne_result$Y[, 1], Y = tsne_result$Y[, 2])

      ggplot(tsne_data, aes(x = X, y = Y, color = as.factor(cluster_allocation))) +
        geom_point(alpha = 0.4) +
        theme_minimal() +
        scale_color_discrete(name = "Cell cluster") +
        labs(title = title)
    }


    png(file.path(output_path, paste0("biclust_data_gen_complete_tsne_",plot_suffix,".png")),
        width = 1024, height = 480, units = "px")
    tsne_plot(indata = dat,
              cluster_allocation = true_cluster_allocation,
              title = "tSNE plot of entire data set ")-> p
    print(p)
    dev.off()

    png(file.path(output_path, paste0("biclust_data_gen_regulators_tsne_",plot_suffix,".png")),
        width = 1024, height = 480, units = "px")
    tsne_plot(indata = Z_r,
              cluster_allocation = true_cluster_allocation,
              title = "tSNE plot of regulators")-> p
    print(p)
    dev.off()

    png(file.path(output_path, paste0("biclust_data_gen_targets_tsne_",plot_suffix,".png")),
        width = 1024, height = 480, units = "px")
    tsne_plot(indata = Z_t,
              cluster_allocation = true_cluster_allocation,
              title = "tSNE plot of targets")-> p
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
    ggplot(df_combined, aes(x = regulator,
                            y = as.factor(target_gene_cluster),
                            fill = as.factor(value))
          )+
      geom_tile() +
      facet_wrap(~ plot, ncol = 2) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            plot.title = element_text(hjust = 0.5)) +
      scale_fill_manual( values = c(
                                    "-1" = "#075AFF",
                                    "0" = "#FFFFCC",
                                    "1" = "#FF0000"
                                    )
                         )+
      xlab("Regulator gene") +
      ylab("Target gene cluster")  +
      ggtitle("Heatmap of S matrices") -> p


    png(file.path(output_path, paste0("biclust_data_gen_S_",plot_suffix,".png")),
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

    ggplot(df_combined, aes(
                            x = target_gene,
                            y = as.factor(target_gene_cluster),
                            fill = as.factor(value)
                            )
          ) +
      geom_tile() +
      facet_wrap(~ plot, ncol = 1) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            plot.title = element_text(hjust = 0.5)) +
      scale_fill_manual(
        values = c("0" = "#075AFF", "1" = "#FF0000")
                       ) +
      xlab("Target gene") +
      ylab("Target gene cluster")  +
      ggtitle("Heatmap of Pi matrices") -> p

    png(file.path(output_path, paste0("biclust_data_gen_Pi_",plot_suffix,".png")),
          width = 1024, height = 480, units = "px")
    print(p)
    dev.off()

    ## #use hclust to find good way to order things
    # dd <- Reduce('+', X)
    # hh <- hclust(as.dist(max(dd)-dd))
    #
    # #Retrieve what to plot
    # listOfMatrices<-lapply(X,function(x) x[hh$order,hh$order])

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
                            y = target_gene,
                            fill = value)) +
      geom_tile() +
      facet_wrap(~ plot, ncol = 2) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            plot.title = element_text(hjust = 0.5)) +
      scale_fill_gradient2(low = "#075AFF",
                           mid = "#FFFFCC",
                           high = "#FF0000")-> p

    png(file.path(output_path, paste0("biclust_data_gen_B_",plot_suffix,".png")),
        width = 1024, height = 480, units = "px")
    print(p)
    dev.off()

    cat("DONE plotting stuff for the simulated cell cluestering data")


  }

  true_Pi      <- lapply(dummy_data, function(x) x$Pi)
  true_S <- lapply(dummy_data, function(x) x$S)

  true_Pi_numbers <- vector(mode = "list", length = length(true_Pi))
  for(i in 1:length(true_Pi)){
    true_Pi_numbers[[i]] <- apply(true_Pi[[i]], 2, function(col) which.max(col != 0))
    # new_target_gene_indexes <- order(true_Pi_numbers[[i]])
    # true_Pi_numbers[[i]] <- true_Pi_numbers[[i]][new_target_gene_indexes]
    # n_target_genes <- max(new_target_gene_indexes)
    # dat[true_cluster_allocation==i, 1:n_target_genes] <- dat[true_cluster_allocation==i, new_target_gene_indexes]
  }

  return(list(
              disturbed_initial_cell_clust = disturbed_initial_cell_clust,
              initial_cell_clust = initial_cell_clust,
              true_cell_clust = true_cluster_allocation,
              true_betas = betas,
              dat = dat,
              true_target_gene_allocation = true_Pi_numbers,
              true_Pi = true_Pi,
              true_S = true_S,
              counts = counts

              )
         )
}




# runs only when script is run by itself
# || interactive()
if (sys.nframe() == 0) {
  # ... do main stuff

  # n_target_genes = 50  # Number of target genes
  # n_regulator_genes = 30  # Number of regulator genes
  # n_cells = 10000 # Number of cells
  # n_target_gene_clusters = 3  # Number of target gene clusters
  # regulator_mean = 1  # Mean expression of regulator genes
  # coefficient_mean = c(5, 20, 100)  # Mean coefficients/betas in true model, length n_target_gene_clusters)
  # coefficient_sd = c(0.1, 0.1, 0.1)
  #
  set.seed(1234) #gives R index close to 1
  res <- generate_dummy_data_for_cell_clustering ()
  print(str(res))
}

