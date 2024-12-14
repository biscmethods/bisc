# This plots which regulators are affecting which target gene modules in the generated data
# Basically we want to translate the guessed labels to the true labels in our guesses to start with
# (in a cell cluster we guess which target genes belong to which gene module)
# After that we move the rows in e.g. true_S (or the corresponding matrix in our guess) so that
# the rows are in the same order in the guess and the generated data.
# Then we can plot comparative networks graphs?
# OR we can calculate the accuracy of our guesses by comparing the true labels to the guessed labels
# E.g. by ROC.

# Network plot of generated data, which regulator genes affect each target gene module ----------------------------------------------------------------------------------------------------------------------------------------------------


# Load required libraries
library(igraph)
library(ggplot2)
library(ggraph)
library(ggforce)
library(concaveman)

scenarios_to_plot <- c(1,12)
set.seed(51124)

# Function to create an edge list with color-coded edges
create_colored_edge_list <- function(true_S_list) {
  # Combine all edge lists from multiple cell clusters
  all_edge_lists <- lapply(true_S_list, function(true_S_matrix) {
    # Get non-zero entries
    non_zero_entries <- which(true_S_matrix != 0, arr.ind = TRUE)

    # Create edge data frame
    edge_df <- data.frame(
      from = paste0("Regulator ", non_zero_entries[, 2]),
      to = paste0("Target Gene Module ", non_zero_entries[, 1]),
      weight = true_S_matrix[non_zero_entries],
      color = ifelse(true_S_matrix[non_zero_entries] > 0, "red", "blue"),
      stringsAsFactors = FALSE
    )

    return(edge_df)
  })

  # Combine edge lists if multiple cell clusters
  edge_list <- do.call(rbind, all_edge_lists)

  return(edge_list)
}


# Create a function to process the data and create an edge list
create_edge_list <- function(data) {
  edge_list <- data.frame()
  for (cluster in seq_along(data)) {
    for (module in 1:nrow(data[[cluster]])) {
      for (gene in 1:ncol(data[[cluster]])) {
        if (data[[cluster]][module, gene] != 0) {
          edge_list <- rbind(edge_list, data.frame(
            from = paste0("Regulator ", gene),
            # to = paste0("Module_", module, "_Cluster_", cluster)
            to = paste0("Target Gene Module ", module)
          ))
        }
      }
    }
  }
  return(edge_list)
}


create_group_list <- function(true_S){
  group_list <- vector(mode="list", length=length(true_S))

  for (cluster in seq_along(true_S)) {
    true_S[[cluster]][true_S[[cluster]]==-1] <- 1
    ones <- true_S[[cluster]]==1
    negative_ones <- true_S[[cluster]]==-1
    ones <- which_rows(ones)
    negative_ones <- which_rows(negative_ones)
    for(i_module in seq_along(ones)){

      if(length(ones[[i_module]])>1){
        group_list[[cluster]][paste0("Regulator group ", i_module)] <- list(ones[[i_module]])
      }
    }
    for(i_module in seq_along(negative_ones)){
      if(length(negative_ones[[i_module]])>1){
        group_list[[cluster]][paste0("Regulator group ", i_module)] <- list(negative_ones[[i_module]])
      }
    }

    # for (module in 1:nrow(true_S[[cluster]])) {
    #   for (gene in 1:ncol(true_S[[cluster]])) {
    #     if (data[[cluster]][module, gene] != 0) {
    #       edge_list <- rbind(edge_list, data.frame(
    #         from = paste0("Gene_", gene),
    #         to = paste0("Module_", module, "_Cluster_", cluster)
    #       ))
    #     }
    #   }
    # }
  }
  return(group_list)
}

which_rows <- function(mat) {
  lapply(seq_len(nrow(mat)), function(i) {
    row_log <- if (sum(mat[i,]) == 1) rep(FALSE, ncol(mat)) else mat[i,]
    row_true <- which(mat[i,])

  })
}


for(i_scenario in scenarios_to_plot){
  # true_S is a list with n_cell_cluster matrices. Each matrix is n_module rows and n_regulator cols.
  # The matrices have 1, -1 or 0. E.g. 1 in (1,3) if cell module 3 is affected positively by regulator 1.
  # [[1]]
  # [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14] [,15] [,16] [,17] [,18] [,19] [,20]
  # [1,]    1    1    1    1    1    1    1    1    0     0     0     0     0     0     0     0     0     0     0     0
  # [2,]    1    1    1    1    0    0    0    0    1     1     1     1     1     1     0     0     0     0     0     0
  # [3,]    1    0    0    0    1   -1    0    0    1     0     0     0     0     0     1     1     1     0     0     0
  #
  # [[2]]
  # [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14] [,15] [,16] [,17] [,18] [,19] [,20]
  # [1,]    1    1    1    1    1    1   -1    1    1    -1     0     0     0     0     0     0     0     0     0     0
  # [2,]    1    1    1   -1   -1    1   -1    1    0     0     1     1     1    -1     1     0     0     0     0     0
  # add summed regulator weights for target gene module? see generate_dummy_data_for_scregclust row 400


  for(i_cell_cluster in seq_along(scenarios[[i_scenario]]$generated_data$true_S)){
    pdf(file.path(output_path, paste0("network_scenario_", i_scenario,"_cell_cluster",i_cell_cluster,".pdf")), width =  5.8, height = 5.8)

    # Get current true_S matrix
    true_S <- scenarios[[i_scenario]]$generated_data$true_S[[i_cell_cluster]]

    # Get gene allocation for this cell cluster
    gene_allocation <- table(scenarios[[i_scenario]]$generated_data$true_target_gene_allocation[[i_cell_cluster]])


    # Create edge list with colors
    edge_list <- create_colored_edge_list(list(true_S))

    group_list <- create_group_list(list(true_S))
    group_list <- unlist(group_list, recursive=FALSE)


    # Create the graph
    edge_list_pos <- edge_list  # We create a positive only version because the layouting doesn't work with negative weights
    edge_list_pos$weight <- 1
    g_pos <- igraph::graph_from_data_frame(edge_list_pos, directed = TRUE)
    g <- igraph::graph_from_data_frame(edge_list, directed = TRUE)


    # Set node types (gene or module)
    V(g)$type <- ifelse(grepl("Regulator ", V(g)$name), "Regulator ", "Target Gene Module ")
    V(g_pos)$type <- ifelse(grepl("Regulator ", V(g)$name), "Regulator ", "Target Gene Module ")

    # Scale module node sizes based on number of genes
    module_sizes <- sapply(names(gene_allocation), function(module_num) {
      module_name <- paste0("Target Gene Module ", module_num)
      if(module_name %in% V(g)$name) {
        # Scale between 5 and 20 based on gene count
        5 + (gene_allocation[module_num] / max(gene_allocation)) * 15
      } else {
        0
      }
    })

    names(module_sizes) <- names(gene_allocation)

    # Create vertex sizes vector
    vertex_sizes <- ifelse(V(g)$type == "Regulator ", 10,
                           sapply(V(g)$name, function(name) {
                             if(grepl("Target Gene Module ", name)) {
                               module_num <- as.numeric(sub("Target Gene Module ", "", name))
                               print(module_num)
                               module_sizes[as.character(module_num)]
                             } else {
                               10
                             }
                           }))

    # Prepare group list for geom_mark_rect
    # Convert group names to node names
    group_ids <- lapply(group_list, function(x) V(g)$name[x])
    # group_ids$"Modules" <- V(g)$name[grepl("Module",V(g)$name)]

    group_color <- brewer.pal(length(group_ids), 'Set1')
    # the fill gets an additional alpha value for transparency:
    group_color_fill <- paste0(group_color, '20')

    plot_layout <- igraph::layout_with_fr(g_pos)  # Fruchterman-Reingold layout

    if(length(group_ids)>0){
      par(mar = rep(0.1, 4))   # reduce margins
    }
    plot(g,
         vertex.color = ifelse(V(g)$type == "Regulator ", "lightblue", "lightgreen"),
         vertex.size = vertex_sizes,
         vertex.frame.color = NA,  # Remove vertex border
         edge.color = E(g)$color,
         edge.arrow.size = 0.5,
         mark.groups = group_ids,
         mark.col = group_color_fill,
         mark.border = group_color,
         layout = plot_layout)
    # if(length(group_ids)>0){
    #   legend('topright', legend = names(group_ids),
    #          col = group_color,
    #          pch = 15, bty = "n",  pt.cex = 1.5, cex = 0.8,
    #          text.col = "black", horiz = FALSE)
    # }
    dev.off()
  }
}
