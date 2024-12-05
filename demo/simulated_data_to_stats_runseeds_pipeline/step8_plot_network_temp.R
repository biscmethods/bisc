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

for(i_scenario in scenarios_to_plot){
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

  which_rows <- function(mat) {
    lapply(seq_len(nrow(mat)), function(i) {
      row_log <- if (sum(mat[i,]) == 1) rep(FALSE, ncol(mat)) else mat[i,]
      row_true <- which(mat[i,])

    })
  }

  create_group_list <- function(true_S){
    group_list <- vector(mode="list", length=length(true_S))
    for (cluster in seq_along(true_S)) {
      ones <- true_S[[cluster]]==1
      negative_ones <- true_S[[cluster]]==-1
      ones <- which_rows(ones)
      negative_ones <- which_rows(negative_ones)
      for(i_module in seq_along(ones)){

        if(length(ones[[i_module]])>1){
          group_list[[cluster]][paste0("Positive regulator module ", i_module)] <- list(ones[[i_module]])
        }
      }
      for(i_module in seq_along(negative_ones)){
        if(length(negative_ones[[i_module]])>1){
          group_list[[cluster]][paste0("Negative regulator module ", i_module)] <- list(negative_ones[[i_module]])
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
    pdf(file.path(output_path, paste0("network_scenario_", i_scenario,"_cell_cluster",i_cell_cluster,".pdf")), width =  7, height = 7)

    edge_list <- create_edge_list(list(scenarios[[i_scenario]]$generated_data$true_S[[i_cell_cluster]]))
    group_list <- create_group_list(list(scenarios[[i_scenario]]$generated_data$true_S[[i_cell_cluster]]))
    group_list <- unlist(group_list, recursive=FALSE)


    # Create the graph
    g <- igraph::graph_from_data_frame(edge_list, directed = FALSE)

    # Set node types (gene or module)
    V(g)$type <- ifelse(grepl("Gene", V(g)$name), "Gene", "Module")

    # Prepare group list for geom_mark_rect
    # Convert group names to node names
    group_ids <- lapply(group_list, function(x) V(g)$name[x])
    # group_ids$"Modules" <- V(g)$name[grepl("Module",V(g)$name)]

    group_color <- brewer.pal(length(group_ids), 'Set1')
    # the fill gets an additional alpha value for transparency:
    group_color_fill <- paste0(group_color, '20')

    plot_layout <- igraph::layout_with_fr(g)  # Fruchterman-Reingold layout

    if(length(group_ids)>0){
      par(mar = rep(0.1, 4))   # reduce margins
    }
    plot(g, vertex.color = 'white', vertex.size = 9,
         edge.color = rgb(0.5, 0.5, 0.5, 0.2),
         mark.groups = group_ids,
         mark.col = group_color_fill,
         mark.border = group_color,
         layout = plot_layout)
    if(length(group_ids)>0){
      legend('topright', legend = names(group_ids),
             col = group_color,
             pch = 15, bty = "n",  pt.cex = 1.5, cex = 0.8,
             text.col = "black", horiz = FALSE)
    }
    dev.off()
  }
}
