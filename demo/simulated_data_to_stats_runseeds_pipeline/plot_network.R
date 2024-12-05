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

i_scenario <- 12

# Create a function to process the data and create an edge list
create_edge_list <- function(data) {
  edge_list <- data.frame()
  for (cluster in seq_along(data)) {
    for (module in 1:nrow(data[[cluster]])) {
      for (gene in 1:ncol(data[[cluster]])) {
        if (data[[cluster]][module, gene] != 0) {
          edge_list <- rbind(edge_list, data.frame(
            from = paste0("Gene_", gene),
            to = paste0("Module_", module, "_Cluster_", cluster)
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
edge_list <- create_edge_list(scenarios[[i_scenario]]$generated_data$true_S)
group_list <- create_group_list(scenarios[[i_scenario]]$generated_data$true_S)
group_list <- unlist(group_list, recursive=FALSE)



# true_Pi  which target genes are in which modules.
# It's a list of n_cell_clusters with matrices. Each matrix is n_modules rows and n_target_genes cols. It has a 1
# for which module each target gene belongs to. It can have only one 1 per column for each cell cluster matrix.
# scenarios[[11]]$generated_data$true_Pi
# scenarios[[11]]$generated_data$true_target_gene_allocation is the same as true_Pi but is not binary but just index instead (e.g. 1111223333)

# Create the graph
g <- igraph::graph_from_data_frame(edge_list, directed = FALSE)

# Set node types (gene or module)
V(g)$type <- ifelse(grepl("Gene", V(g)$name), "Gene", "Module")


# Prepare group list for geom_mark_rect
# Convert group names to node names
group_ids <- lapply(group_list, function(x) V(g)$name[x])

group_color <- brewer.pal(length(group_ids), 'Set1')
# the fill gets an additional alpha value for transparency:
group_color_fill <- paste0(group_color, '20')

# Use the same layout coordinates for both the igraph plot and the polygon extraction
plot_layout <- igraph::layout_with_fr(g)
plot_layout <- igraph::norm_coords(plot_layout, ymin=0, ymax=1, xmin=0, xmax=1)
plot_layout[21,] = c(0.3, 0.3)
#
# # Plot the igraph graph
# Create the plot
p <- ggraph(g, layout = plot_layout) +
  geom_edge_link(alpha = 0.4) +
  geom_node_point(aes(color = type, size = type)) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_color_manual(values = c("Gene" = "lightblue", "Module" = "lightgreen")) +
  scale_size_manual(values = c("Gene" = 5, "Module" = 7)) +
  theme_minimal() +
  theme(legend.position = "bottom")

# print(p)

extract_mark_group_polygons <- function(g, group_ids, plot_layout) {
  # Normalize the layout coordinates
  x_min <- min(plot_layout[, 1])
  x_max <- max(plot_layout[, 1])
  y_min <- min(plot_layout[, 2])
  y_max <- max(plot_layout[, 2])

  plot_layout_norm <- cbind((plot_layout[, 1] - x_min) / (x_max - x_min),
                            (plot_layout[, 2] - y_min) / (y_max - y_min))

  polygon_coords <- lapply(group_ids, function(group) {
    # Get the indices of nodes in the group
    group_nodes <- which(V(g)$name %in% group)

    # Extract the normalized layout coordinates for these nodes
    node_coords <- plot_layout_norm[group_nodes, ]

    # Use convex hull to create a smooth border polygon
    hull_indices <- chull(node_coords)
    hull_coords <- node_coords[hull_indices, ]

    return(hull_coords)
  })

  return(polygon_coords)
}

# Extract polygon coordinates
mark_group_polygons <- extract_mark_group_polygons(g, group_ids, plot_layout)



# Function to add mark group polygons to the ggplot2 plot
add_mark_group_polygons <- function(p, mark_group_polygons) {
  for (i in seq_along(mark_group_polygons)) {
    polygon_coords <- as.data.frame(mark_group_polygons[[i]])
    colnames(polygon_coords) <- c("x", "y")
    p <- p + geom_polygon(data = polygon_coords,
                          aes(x = x, y = y),
                          fill = group_color_fill[i],
                          color = group_color[i],
                          alpha = 0.2)
  }
  return(p)
}

# Add the mark group polygons to the ggplot2 plot
p <- add_mark_group_polygons(p, mark_group_polygons)

# Print the final plot
print(p)
